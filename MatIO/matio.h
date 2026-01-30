// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2025, the Anboto author and contributors
#ifndef _plugin_matio_matio_h_
#define _plugin_matio_matio_h_

#include <plugin/matio/lib/matio.h>
#include <plugin/matio/lib/matioConfig.h>
#include <Eigen/Eigen.h>
#include <Eigen/MultiDimMatrix.h>


namespace Upp {

struct MatioError : Exc {
	MatioError(const String& desc) : Exc("MatIO error: " + desc) {}
};

inline void Verify(bool cond, const char* msg = "Problem found") { if (!cond) throw MatioError(msg); }
inline void VerifyError(int err, String msg = "Problem found") { 
	if (err != 0) {
		if (err == MATIO_E_GENERIC_WRITE_ERROR)
			msg << ". " << "The .mat file may be already opened";
		throw MatioError(msg); 
	}
}


class MatVar {
public:
	MatVar() : var(nullptr) {}
	MatVar(matvar_t *var_, bool del_ = true) : var(var_), del(del_) {}
	~MatVar() {
		if (var && del)
			Mat_VarFree(var);
		var = nullptr;
	}
	operator matvar_t*() { return var; }
    
	const char *GetName()			{Verify(var, "Variable is not created"); return var->name;}
	enum matio_classes GetType() 	{Verify(var, "Variable is not created"); return var->class_type;}
	int GetDimCount() const			{Verify(var, "Variable is not created"); return var->rank;}
	int GetDimCount(int dim) const	{Verify(var, "Variable is not created"); return (int)var->dims[dim];}
	
	const char* GetTypeString()	{
		Verify(var, "Variable is not created"); 
		
	    const char *class_type_desc[] = {"Undefined","Cell Array","Structure",
	       "Object","Character Array","Sparse Array","Double Precision Array",
	       "Single Precision Array", "8-bit, signed integer array",
	       "8-bit, unsigned integer array","16-bit, signed integer array",
	       "16-bit, unsigned integer array","32-bit, signed integer array",
	       "32-bit, unsigned integer array","64-bit, signed integer array",
	       "64-bit, unsigned integer array","Function","Opaque"};
	
		return class_type_desc[var->class_type];
	}
	
	int GetCount() {
		Verify(var, "Variable is not created"); 
		
		int ret = GetDimCount(0);
		for (int i = 1; i < GetDimCount(); ++i)
			ret *= GetDimCount(i);
		return ret;
	}
	
	int GetFieldCount() {
		Verify(var, "Variable is not created"); 
		Verify(IsStruct(), "Variable is not a struct");
		
		return (int)Mat_VarGetNumberOfFields(var);
	}
	const char *GetFieldName(int id) {	
		Verify(var, "Variable is not created"); 
		Verify(IsStruct(), "Variable is not a struct");
		
		return Mat_VarGetStructFieldnames(var)[id];
	}

	bool IsLoaded() {return var != 0;}
	bool IsStruct() {return var->class_type == MAT_C_STRUCT;}
	
	MatVar GetStructVar(const char *field_name) {
		return MatVar(Mat_VarGetStructFieldByName(var, field_name, 0), false);	
	}
	bool ExistStructVar(const char *field_name) {
		return Mat_VarGetStructFieldByName(var, field_name, 0);	
	}
	
protected:
	matvar_t *var = nullptr;
	bool del = true;
};

static inline bool IsScalar(const matvar_t* v) {
    if (!v) 
    	return false;
    size_t n = 1; 
    for (int i = 0; i < v->rank; ++i) 
    	n *= v->dims[i];
    return n == 1;
}


template <class T>
static inline Vector<T> RowToCol(const T* row, size_t rows, size_t cols) {
    Vector<T> col((int)(rows*cols));
    for (size_t r = 0; r < rows; ++r)
        for (size_t c = 0; c < cols; ++c)
            col[(int)(c*rows + r)] = row[(int)(r*cols + c)];
    return col;
}

template <class T>
static inline Vector<T> ColToRow(const T* col, size_t rows, size_t cols) {
    Vector<T> row((int)(rows*cols));
    for (size_t c = 0; c < cols; ++c)
        for (size_t r = 0; r < rows; ++r)
            row[(int)(r*cols + c)] = col[(int)(c*rows + r)];
    return row;
}


class MatFile {
public:
    MatFile()  {}
    ~MatFile() {Close(); }

    MatFile(const MatFile&) = delete;
    MatFile& operator=(const MatFile&) = delete;

    MatFile(MatFile&& o) noexcept {Move(pick(o)); }
    MatFile& operator=(MatFile&& o) noexcept { 
    	if (this != &o) {
    		Close(); 
    		Move(pick(o)); 
    	} 
    	return *this; 
    }

    bool OpenRead(const char *path) {
        path_ = path;
        fp_ = Mat_Open(path, MAT_ACC_RDONLY);
        return fp_;
    }
    bool OpenReadWrite(const char *path) {
        path_ = path;
        fp_ = Mat_Open(path, MAT_ACC_RDWR);
        return fp_;
    }
    bool OpenCreate(const char *path, mat_ft ver = mat_ft::MAT_FT_MAT5) {
        path_ = path;
        fp_ = Mat_CreateVer(path, nullptr, ver);
        return fp_;
    }

	mat_ft GetVersion() {
		Verify(fp_, "File is not created");
		
		return Mat_GetVersion(fp_);
	}
	String GetVersionName() {
		mat_ft ver = GetVersion();
		switch (ver) {
		case MAT_FT_MAT4: 		return "4";
		case MAT_FT_MAT5: 		return "5";
		case MAT_FT_MAT73:		return "7.3";
		case MAT_FT_UNDEFINED:	return "unknown";
		}
		return "unknown";
	}
		
	Vector<String> GetVarList() const {
		Verify(fp_, "File is not created");
		
		size_t numVar = 0; 
		char *const* list = Mat_GetDir(fp_, &numVar);
		
		Vector<String>ret((int)numVar);
		for (size_t i = 0; i < numVar; ++i)
			ret[(int)i] = list[i];
		
		return ret;
	}
	
	String GetVarName(const char* name) const {
		String namelower = ToLower(name);
		Vector<String> list = GetVarList();
		for (const String &namereal : list)
			if (ToLower(namereal) == namelower)
				return namereal;
		return String();
	}
	
    mat_t* Mat_t() const noexcept { return fp_; }

    
    MatVar GetVar(const char* name, bool nocase = false) const {
        if (nocase) {
            String rname = GetVarName(name);
            return GetVar(rname);
        }
        Verify(fp_, "File not open");
        matvar_t* v = Mat_VarRead(fp_, name);
        if (v)
        	return MatVar(v);
        return MatVar();
    }
    MatVar GetVar(const Vector<String> &names, bool nocase = false) const {
        Verify(fp_, "File not open");
        for (String name : names) {
            if (nocase)
	            name = GetVarName(name);
        	matvar_t* v = Mat_VarRead(fp_, name);
        	if (v)
        		return MatVar(v);
        }
        return MatVar();
    }
    bool Exist(const char* name) const {
        Verify(fp_, "File not open");
        matvar_t* v = Mat_VarRead(fp_, name);
        if (v) {
            Mat_VarFree(v);
            return true;
        } else
            return false;
    }

    
    void Write(const char* name, int v, matio_compression compression = MAT_COMPRESSION_NONE) {
        size_t dims[2] = {1,1};
        int tmp = v; // stack
        matvar_t* mv = Mat_VarCreate(name, MAT_C_INT64, MAT_T_INT64, 2, dims, &tmp, 0);
        Verify(mv, "Write(int): Mat_VarCreate failed");
        VerifyError(Mat_VarWrite(fp_, mv, compression), "Write(int): Mat_VarWrite failed");
        Mat_VarFree(mv);
    }
    void Write(const char* name, double v, matio_compression compression = MAT_COMPRESSION_NONE) {
        size_t dims[2] = {1,1};
        double tmp = v;
        matvar_t* mv = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &tmp, 0);
        Verify(mv, "Write(double): Mat_VarCreate failed");
        VerifyError(Mat_VarWrite(fp_, mv, compression), "Write(double): Mat_VarWrite failed");
        Mat_VarFree(mv);
    }
    void Write(const char* name, std::complex<double> z, matio_compression compression = MAT_COMPRESSION_NONE) {
        size_t dims[2] = {1,1};
        double re = z.real(), im = z.imag();       // stack buffers
        mat_complex_split_t split{&re, &im};
        matvar_t* mv = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &split, MAT_F_COMPLEX);
        Verify(mv, "Write(complex): Mat_VarCreate failed");
        VerifyError(Mat_VarWrite(fp_, mv, compression), "Write(complex): Mat_VarWrite failed");
        Mat_VarFree(mv);
    }

    void Write(const char* name, const char *s, matio_compression compression = MAT_COMPRESSION_NONE) {
        String utf8(s);
        size_t dims[2] = {1, (size_t)utf8.GetLength()};
        Vector<uint8_t> bytes((int)dims[1]);
        for (size_t i = 0; i < dims[1]; ++i) 
        	bytes[(int)i] = (uint8_t)utf8[(int)i];
        matio_types tchar =
        #ifdef MAT_T_UTF8
            MAT_T_UTF8;
        #else
            MAT_T_UINT8;
        #endif
        matvar_t* mv = Mat_VarCreate(name, MAT_C_CHAR, tchar, 2, dims, bytes.empty() ? NULL : bytes.begin(), 0);
        Verify(mv, "Write: Mat_VarCreate failed");
        VerifyError(Mat_VarWrite(fp_, mv, compression), "Write: Mat_VarWrite failed");
        Mat_VarFree(mv);
    }

	template <class T>
	void Write(const char* name, const Vector<T> &data, matio_compression compression = MAT_COMPRESSION_NONE) {
        if constexpr (std::is_same<T,double>::value)
            WriteColMajor(name, 1, (size_t)data.size(), data, compression);
        else if constexpr (std::is_same<T,int>::value)
            WriteColMajor(name, 1, (size_t)data.size(), data, compression);
        else
            static_assert(sizeof(T)==0, "Supported T: double or int");
    }                            
	void Write(const char* name, const Vector<std::complex<double>> &data, matio_compression compression = MAT_COMPRESSION_NONE) {
        WriteColMajor(name, 1, (size_t)data.size(), data, compression);
    }
                                     
    
    void WriteColMajor(const char* name, size_t rows, size_t cols,
                             const double* colMajor, matio_compression compression = MAT_COMPRESSION_NONE) {
        size_t dims[2] = {rows, cols};
        Vector<double> buf; 	Copy(colMajor, rows*cols, buf);
        matvar_t* mv = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, buf.Begin(), 0);
        Verify(mv, "WriteColMajor(double): Mat_VarCreate failed");
        VerifyError(Mat_VarWrite(fp_, mv, compression), "WriteMatrixColMajor(double): Mat_VarWrite failed");
        Mat_VarFree(mv);
    }
    void WriteColMajor(const char* name, size_t rows, size_t cols,
                             const int* colMajor, matio_compression compression = MAT_COMPRESSION_NONE) {
        size_t dims[2] = {rows, cols};
        Vector<int> buf;	Copy(colMajor, rows*cols, buf);
        matvar_t* mv = Mat_VarCreate(name, MAT_C_INT64, MAT_T_INT64, 2, dims, buf.Begin(), 0);
        Verify(mv, "WriteColMajor(int): Mat_VarCreate failed");
        VerifyError(Mat_VarWrite(fp_, mv, compression), "WriteMatrixColMajor(int): Mat_VarWrite failed");
        Mat_VarFree(mv);
    }
    template <class T>
    void WriteRowMajor(const char* name, size_t rows, size_t cols,
                             const T* rowMajor, matio_compression compression = MAT_COMPRESSION_NONE) {
        auto col = RowToCol(rowMajor, rows, cols);
        if constexpr (std::is_same<T,double>::value)
            WriteColMajor(name, rows, cols, col, compression);
        else if constexpr (std::is_same<T,int>::value)
            WriteColMajor(name, rows, cols, col, compression);
        else
            static_assert(sizeof(T)==0, "Supported T: double or int");
    }

    
    void WriteColMajor(const char* name, size_t rows, size_t cols, const std::complex<double>* data, matio_compression compression = MAT_COMPRESSION_NONE) {
        size_t dims[2] = {rows, cols};
        size_t N = rows*cols;
  		Vector<double> re((int)N), im((int)N);
  		for (size_t i = 0; i < N; ++i) {
  			re[(int)i] = data[(int)i].real();
  			im[(int)i] = data[(int)i].imag();
  		}
        mat_complex_split_t split{re.Begin(), im.Begin()};
        matvar_t* mv = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &split, MAT_F_COMPLEX);
        Verify(mv, "WriteComplexMatrixColMajor: Mat_VarCreate failed");
        VerifyError(Mat_VarWrite(fp_, mv, compression), "WriteComplexMatrixColMajor: Mat_VarWrite failed");
        Mat_VarFree(mv);
    }
    void WriteRowMajor(const char* name, size_t rows, size_t cols,
                                    const std::complex<double>* data,
                                    matio_compression compression = MAT_COMPRESSION_NONE) {
        auto reCol = RowToCol(data, rows, cols);
        WriteColMajor(name, rows, cols, reCol, compression);
    }

    template <class T>
    void WriteRowMajor(const char* name, size_t rows, size_t cols,
                                const Upp::Vector<T>& rowFlat, matio_compression compression = MAT_COMPRESSION_NONE) {
        Verify((size_t)rowFlat.GetCount() == rows*cols, "WriteMatrixUppRowMajor: count mismatch");
        if constexpr (std::is_same<T,double>::value)
            WriteMatrixRowMajor(name, rows, cols, rowFlat.Begin(), compression);
        else if constexpr (std::is_same<T,int>::value)
            WriteMatrixRowMajor(name, rows, cols, rowFlat.Begin(), compression);
        else
            Verify(sizeof(T)==0, "Supported U++ T: double or int");
    }

    template <class Derived>
    void Write(const char* name, const Eigen::MatrixBase<Derived>& M,
                    matio_compression compression = MAT_COMPRESSION_NONE) {
        using Scalar = typename Derived::Scalar;
        const size_t rows = (size_t)M.rows();
        const size_t cols = (size_t)M.cols();

        if constexpr (std::is_same<Scalar,double>::value) {
            if (Derived::IsRowMajor)
                WriteRowMajor(name, rows, cols, M.derived().data(), compression);
            else
                WriteColMajor(name, rows, cols, M.derived().data(), compression);
        } else if constexpr (std::is_same<Scalar,int>::value) {
            if (Derived::IsRowMajor)
                WriteRowMajor(name, rows, cols, M.derived().data(), compression);
            else
                WriteColMajor(name, rows, cols, M.derived().data(), compression);
        } else if constexpr (std::is_same<Scalar,std::complex<double>>::value) {
            Vector<std::complex<double>> data(rows*cols);
            if (Derived::IsRowMajor) {	// Convert row-major to column-major order
                for (size_t r=0; r<rows; ++r)
                    for (size_t c=0; c<cols; ++c)
                        data[c*rows + r] = M(Eigen::Index(r*cols + c));
            } else {
                for (size_t i=0; i<rows*cols; ++i)
                    data[i] = M.derived()(Eigen::Index(i));
            }
            WriteColMajor(name, rows, cols, data, compression);
        } else
            Verify(sizeof(Scalar)==0, "Unsupported Eigen scalar (use double, int, complex<double>)");
    }

	template <class T>
	void Write(const char* name, const MultiDimMatrix<T>& M, matio_compression compression = MAT_COMPRESSION_NONE) {
	    const Vector<int>& ad = M.GetAxisDim();
	    Verify(!ad.IsEmpty(), "Write(MultiDimMatrix): empty dims");            
	    Vector<size_t> dims = ToSizeDims(ad);
	    const size_t N = (size_t)M.size();
	    Verify(N == (size_t)Product(dims), "Write(MultiDimMatrix): size mismatch");
	
	    if constexpr (std::is_same<T, double>::value) {
	        Vector<double> buf; buf.SetCount((int)N);
	        std::memcpy(buf.Begin(), M.begin(), N * sizeof(double));
	        matvar_t* mv = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE,
	                                     dims.size(), dims.Begin(), buf.Begin(), 0);
	        Verify(mv, "Write(MultiDimMatrix<double>): Mat_VarCreate failed");
	        VerifyError(Mat_VarWrite(fp_, mv, compression), "Write(MultiDimMatrix<double>): Mat_VarWrite failed");
	        Mat_VarFree(mv);
	    } else if constexpr (std::is_same<T, int>::value) {
	        Vector<int> buf; buf.SetCount((int)N);
	        std::memcpy(buf.Begin(), M.begin(), N * sizeof(int));
	        matvar_t* mv = Mat_VarCreate(name, MAT_C_INT64, MAT_T_INT64,
	                                     dims.size(), dims.Begin(), buf.Begin(), 0);
	        Verify(mv, "Write(MultiDimMatrix<int>): Mat_VarCreate failed");
	        VerifyError(Mat_VarWrite(fp_, mv, compression), "Write(MultiDimMatrix<int>): Mat_VarWrite failed");
	        Mat_VarFree(mv);
	    } else if constexpr (std::is_same<T, std::complex<double>>::value) {
	        Vector<double> re, im;
	        re.SetCount((int)N); im.SetCount((int)N);
	        const std::complex<double>* p = M.begin();
	        for (size_t i = 0; i < N; ++i) { re[(int)i] = p[i].real(); im[(int)i] = p[i].imag(); }
	        mat_complex_split_t split{re.Begin(), im.Begin()};
	        matvar_t* mv = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE,
	                                     dims.size(), dims.Begin(), &split, MAT_F_COMPLEX);
	        Verify(mv, "Write(MultiDimMatrix<complex>): Mat_VarCreate failed");
	        VerifyError(Mat_VarWrite(fp_, mv, compression), "Write(MultiDimMatrix<complex>): Mat_VarWrite failed");
	        Mat_VarFree(mv);
	    } else
	        static_assert(sizeof(T) == 0, "Unsupported MultiDimMatrix<T> (use double, int, std::complex<double>)");
	}

	void WriteCell(const char* name, const Vector<String> &strs, matio_compression compression = MAT_COMPRESSION_NONE) {
		size_t num = (size_t)strs.size();
		Buffer<matvar_t *> cell_elements(num);
	
	    for (int i = 0; i < (int)num; i++) {
	        const char *c_str = strs[i];
	
	        size_t dims[2] = {1, strlen(c_str)};
	
	        cell_elements[i] = Mat_VarCreate(NULL, MAT_C_CHAR, MAT_T_UTF8, 2, dims, c_str, MAT_F_DONT_COPY_DATA);
	    }
	    size_t cell_dims[2] = {num, 1};
	
	    matvar_t *cell_array = Mat_VarCreate(name, MAT_C_CELL, MAT_T_CELL, 2, cell_dims, cell_elements, MAT_F_DONT_COPY_DATA);
	
	    Mat_VarWrite(fp_, cell_array, compression);
	
	    for (size_t i = 0; i < num; i++) 
	        Mat_VarFree(cell_elements[i]);
	    
	    Mat_VarFree(cell_array);
	}
	
	class StructNode {
	public:
	    StructNode() {}
	    StructNode(MatFile &mf, const char* name, const Vector<String>& fieldNames) : mf_(&mf) {Create(name, fieldNames);}
	    void Create(const char* name, const Vector<String>& fieldNames) {
	        name_ = name;
	        Vector<const char*> fields; fields.Reserve(fieldNames.size());
	        for (const auto& s : fieldNames) 
	        	fields << s;
	        size_t dims[2] = {1,1};
	        var_ = Mat_VarCreateStruct(name, 2, dims, fields, (unsigned)fields.size());
	        Verify(var_, "StructNode: Mat_VarCreateStruct failed");
	    }
	
	    StructNode(const StructNode&) = delete;
	    StructNode& operator=(const StructNode&) = delete;
	
	    StructNode(StructNode&& other) noexcept { Move(pick(other)); }
	    StructNode& operator=(StructNode&& other) noexcept {
	        if (this != &other) { 
	        	Clear(); 
	        	Move(pick(other)); 
	        }
	        return *this;
	    }
	    ~StructNode() {Clear();}
	
	    void Write(const char* field, int v) {
	        size_t dims[2] = {1,1};
	        scalars_i64_ << v;
	        int* p = &scalars_i64_.back();
	        matvar_t* child = Mat_VarCreate(field, MAT_C_INT64, MAT_T_INT64, 2, dims, p, 0);
	        Verify(child, "StructNode Set(int): Mat_VarCreate failed");
	        Verify(Mat_VarSetStructFieldByName(var_, field, 0, child) == 0,
	               "StructNode Set(int): set field failed");
	    }
	
	    void Write(const char* field, double v) {
	        size_t dims[2] = {1,1};
	        scalars_double_ << v;
	        double* p = &scalars_double_.back();
	        matvar_t* child = Mat_VarCreate(field, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, p, 0);
	        Verify(child, "StructNode Set(double): Mat_VarCreate failed");
	        Verify(Mat_VarSetStructFieldByName(var_, field, 0, child) == 0,
	               "StructNode Set(double): set field failed");
	    }
	
	    void Write(const char* field, const String& s) {
	        const int N = s.GetLength();
	        Vector<uint8_t> &nw = strings_.Add();
	        nw.SetCount(N);
	        for (int i=0;i<N;++i) 
	        	nw[i] = (uint8_t)s[i];
	        size_t dims[2] = {1, (size_t)N};
	        matio_types tchar =
	        #ifdef MAT_T_UTF8
	            MAT_T_UTF8;
	        #else
	            MAT_T_UINT8;
	        #endif
	        matvar_t* child = Mat_VarCreate(field, MAT_C_CHAR, tchar, 2, dims, strings_.back(), 0);
	        Verify(child, "StructNode Set(String): Mat_VarCreate failed");
	        Verify(Mat_VarSetStructFieldByName(var_, field, 0, child) == 0,
	               "StructNode Set(String): set field failed");
	    }
	
	    void Write(const char* field, std::complex<double> z) {
	        complex_re << z.real();
	        complex_im << z.imag();
	        
	        size_t dims[2] = {1,1};
	        mat_complex_split_t split{&complex_re.back(), &complex_im.back()};
	        matvar_t* child = Mat_VarCreate(field, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &split, MAT_F_COMPLEX);
	        Verify(child, "StructNode Set(complex): Mat_VarCreate failed");
	        Verify(Mat_VarSetStructFieldByName(var_, field, 0, child) == 0,
	               "StructNode Set(complex): set field failed");
	    }
	
	    void WriteColMajor(const char* field, size_t rows, size_t cols, const double* colMajor) {
	        Vector<double> &nw = mats_double_.Add();
	        Copy(colMajor, rows*cols, nw);
	        size_t dims[2] = {rows, cols};
	        matvar_t* child = Mat_VarCreate(field, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, mats_double_.back(), 0);
	        Verify(child, "StructNode SetM: Mat_VarCreate failed");
	        Verify(Mat_VarSetStructFieldByName(var_, field, 0, child) == 0,
	               "StructNode Set: set field failed");
	    }
	
	    void WriteColMajor(const char* field, size_t rows, size_t cols, const int* colMajor) {
	        Vector<int> &nw = mats_i64_.Add();
	        Copy(colMajor, rows*cols, nw);
	        size_t dims[2] = {rows, cols};
	        matvar_t* child = Mat_VarCreate(field, MAT_C_INT64, MAT_T_INT64, 2, dims, mats_i64_.back(), 0);
	        Verify(child, "StructNode SetMatrixInt: Mat_VarCreate failed");
	        Verify(Mat_VarSetStructFieldByName(var_, field, 0, child) == 0,
	               "StructNode SetMatrixInt: set field failed");
	    }
	
	    void WriteColMajor(const char* field, size_t rows, size_t cols, const std::complex<double>* data) {
	        Vector<double> &nwre = mats_cplx_re_.Add(),
	        			   &nwim = mats_cplx_im_.Add();
	        int N = (int)(rows*cols);
	        nwre.SetCount(N);
	        nwim.SetCount(N);
	        for (int i = 0; i < N; ++i) {
	            nwre[i] = data[i].real();
	            nwim[i] = data[i].imag();
	        }
	        mat_complex_split_t split{mats_cplx_re_.back(), mats_cplx_im_.back()};
	        size_t dims[2] = {rows, cols};
	        matvar_t* child = Mat_VarCreate(field, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &split, MAT_F_COMPLEX);
	        Verify(child, "StructNode SetComplexMatrix: Mat_VarCreate failed");
	        Verify(Mat_VarSetStructFieldByName(var_, field, 0, child) == 0,
	               "StructNode SetComplexMatrix: set field failed");
	    }
	    
	    template <class T>
		void Write(const char* name, const Vector<T> &data) {
	        if constexpr (std::is_same<T,double>::value)
	            WriteColMajor(name, 1, (size_t)data.size(), data);
	        else if constexpr (std::is_same<T,int>::value)
	            WriteColMajor(name, 1, (size_t)data.size(), data);
	        else
	            static_assert(sizeof(T)==0, "Supported T: double or int");
	    }                            
		void Write(const char* name, const Vector<std::complex<double>> &data) {
	        WriteColMajor(name, 1, (size_t)data.size(), data);
	    }
	    
		template <class T>
		void Write(const char* field, const MultiDimMatrix<T>& M) {
		    const Vector<int>& ad = M.GetAxisDim();
		    Verify(!ad.IsEmpty(), "StructNode Write(MultiDimMatrix): empty dims");
		    Vector<size_t> dims = ToSizeDims(ad);
		    const size_t N = (size_t)M.size();
		
		    if constexpr (std::is_same<T, double>::value) {
		        Vector<double>& nw = mats_double_.Add();
		        nw.SetCount((int)N);
		        std::memcpy(nw.Begin(), M.begin(), N * sizeof(double));
		        matvar_t* child = Mat_VarCreate(field, MAT_C_DOUBLE, MAT_T_DOUBLE,
		                                        dims.size(), dims.Begin(), mats_double_.back(), 0);
		        Verify(child, "StructNode Write(MultiDimMatrix<double>): Mat_VarCreate failed");
		        Verify(Mat_VarSetStructFieldByName(var_, field, 0, child) == 0,
		               "StructNode Write(MultiDimMatrix<double>): set field failed");
		    } else if constexpr (std::is_same<T, int>::value) {
		        Vector<int>& nw = mats_i64_.Add();
		        nw.SetCount((int)N);
		        std::memcpy(nw.Begin(), M.begin(), N * sizeof(int));
		        matvar_t* child = Mat_VarCreate(field, MAT_C_INT64, MAT_T_INT64,
		                                        dims.size(), dims.Begin(), mats_i64_.back(), 0);
		        Verify(child, "StructNode Write(MultiDimMatrix<int>): Mat_VarCreate failed");
		        Verify(Mat_VarSetStructFieldByName(var_, field, 0, child) == 0,
		               "StructNode Write(MultiDimMatrix<int>): set field failed");
		    } else if constexpr (std::is_same<T, std::complex<double>>::value) {
		        Vector<double>& re = mats_cplx_re_.Add();
		        Vector<double>& im = mats_cplx_im_.Add();
		        re.SetCount((int)N); im.SetCount((int)N);
		        const std::complex<double>* p = M.begin();
		        for (size_t i = 0; i < N; ++i) { re[(int)i] = p[i].real(); im[(int)i] = p[i].imag(); }
		        mat_complex_split_t split{mats_cplx_re_.back(), mats_cplx_im_.back()};
		        matvar_t* child = Mat_VarCreate(field, MAT_C_DOUBLE, MAT_T_DOUBLE,
		                                        dims.size(), dims.Begin(), &split, MAT_F_COMPLEX);
		        Verify(child, "StructNode Write(MultiDimMatrix<complex>): Mat_VarCreate failed");
		        Verify(Mat_VarSetStructFieldByName(var_, field, 0, child) == 0,
		               "StructNode Write(MultiDimMatrix<complex>): set field failed");
		    } else {
		        static_assert(sizeof(T) == 0, "Unsupported MultiDimMatrix<T> (use double, int, std::complex<double>)");
		    }
		}

	
	    StructNode& AddChild(const char* fieldName,
	                             const Vector<String>& childFields)
	    {
	        child_field_names_ << fieldName;
	        StructNode &nw = children_.Add();
	        nw.Create(fieldName, childFields);   // child struct named as the field (MAT accepts it)
	        return children_.back();
	    }
	
	    void Write(matio_compression compression = MAT_COMPRESSION_NONE) {
	        Verify(var_ != nullptr, "StructNode: already written or moved");
	        FinalizeTreeAttachChildren();
	        VerifyError(Mat_VarWrite(mf_->Mat_t(), var_, compression), "StructNode: Mat_VarWrite failed");
	        Mat_VarFree(var_);
	        var_ = nullptr;
	    }
	
	private:
	    void Clear() noexcept { 
	    	if (var_) { 
	    		Mat_VarFree(var_); 
	    		var_ = nullptr; 
	    	} 
	    }
	    void Move(StructNode&& o) noexcept {
	        name_ = pick(o.name_);
	        var_  = o.var_;  o.var_ = nullptr;
	
	        scalars_i64_     = pick(o.scalars_i64_);
	        scalars_double_  = pick(o.scalars_double_);
	        strings_         = pick(o.strings_);
	        complex_re 		 = pick(o.complex_re);
	        complex_im 		 = pick(o.complex_im);
	        mats_double_     = pick(o.mats_double_);
	        mats_i64_        = pick(o.mats_i64_);
	        mats_cplx_re_    = pick(o.mats_cplx_re_);
	        mats_cplx_im_    = pick(o.mats_cplx_im_);
	        children_        = pick(o.children_);
	        child_field_names_ = pick(o.child_field_names_);
	    }
		
	    void FinalizeTreeAttachChildren() {
	        for (int i = 0; i < children_.size(); ++i) {
	            StructNode& ch = children_[i];
	            ch.FinalizeTreeAttachChildren();
	            
	            Verify(Mat_VarSetStructFieldByName(var_, child_field_names_[i], 0, ch.var_) == 0,
	                   "StructNode: attach child failed");
	            ch.var_ = nullptr; // parent now owns child's matvar
	        }
	    }
	    
		String name_;
	    
	    matvar_t* var_ = nullptr;
	    MatFile *mf_ = nullptr;
	
	    Vector<int>             scalars_i64_;
	    Vector<double>          scalars_double_;
	    Vector<Vector<uint8_t>>	strings_;
	    Vector<double>        	complex_re;
	    Vector<double>        	complex_im;
	
	    Vector<Vector<double>>  mats_double_;
	    Vector<Vector<int>>     mats_i64_;
	    Vector<Vector<double>>  mats_cplx_re_;
	    Vector<Vector<double>>  mats_cplx_im_;
	
	    Array<StructNode>       children_;
	    Vector<String>          child_field_names_;
	};

    static int ReadInt(const matvar_t* v) {
        Verify(v && IsScalar(v), "ReadScalarInt: expected scalar");
        switch (v->class_type) {
            case MAT_C_INT64:  return *(int*)v->data;
            case MAT_C_UINT64: {
                unsigned u = *(unsigned*)v->data;
                if (u > (unsigned)std::numeric_limits<int>::max()) 
                	throw MatioError("uint64 exceeds int64");
                return (int)u;
            }
            case MAT_C_DOUBLE: return (int)*(double*)v->data;
            case MAT_C_SINGLE: return (int)*(float*)v->data;
            case MAT_C_INT8:   return (int)*(int8_t*)v->data;
            case MAT_C_UINT8:  return (int)*(uint8_t*)v->data;
            case MAT_C_INT16:  return (int)*(int16_t*)v->data;
            case MAT_C_UINT16: return (int)*(uint16_t*)v->data;
            case MAT_C_INT32:  return (int)*(int32_t*)v->data;
            case MAT_C_UINT32: return (int)*(uint32_t*)v->data;
            default: throw MatioError("ReadScalarInt64: not numeric");
        }
    }
    static double ReadDouble(const matvar_t* v) {
        Verify(v && IsScalar(v), "ReadScalarDouble: expected scalar");
        switch (v->class_type) {
            case MAT_C_DOUBLE: return *(double*)v->data;
            case MAT_C_SINGLE: return (double)*(float*)v->data;
            case MAT_C_INT64:  return (double)*(int64_t*)v->data;
            case MAT_C_UINT64: return (double)*(uint64_t*)v->data;
            case MAT_C_INT32:  return (double)*(int32_t*)v->data;
            case MAT_C_UINT32: return (double)*(uint32_t*)v->data;
            case MAT_C_INT16:  return (double)*(int16_t*)v->data;
            case MAT_C_UINT16: return (double)*(uint16_t*)v->data;
            case MAT_C_INT8:   return (double)*(int8_t*)v->data;
            case MAT_C_UINT8:  return (double)*(uint8_t*)v->data;
            default: throw MatioError("ReadScalarDouble: not numeric");
        }
    }
    static std::complex<double> ReadComplex(const matvar_t* v) {
        Verify(v && IsScalar(v), "ReadScalarComplexDouble: expected scalar");
        Verify(v->class_type == MAT_C_DOUBLE && v->isComplex, "ReadScalarComplexDouble: not complex<double>");
        const auto* cs = static_cast<const mat_complex_split_t*>(v->data);
        const double re = cs ? *(const double*)cs->Re : 0.0;
        const double im = cs ? *(const double*)cs->Im : 0.0;
        return {re, im};
    }

	static String ReadString(const matvar_t* v) {
	    Verify(v && v->class_type == MAT_C_CHAR, "Expected char array");
	    Verify(v->rank == 2, "Expected rank-2");
	    size_t rows = v->dims[0], cols = v->dims[1];
	    size_t N = (rows == 1 ? cols : (cols == 1 ? rows : 0));
	    Verify(N > 0, "Expected 1xN or Nx1 char vector");
	    String s;
	
	    if (v->data_type == MAT_T_UTF8) {
	        const uint8_t* p = static_cast<const uint8_t*>(v->data);
	        for (size_t i=0; i<N; ++i) 
	        	s.Cat((char)p[i]);
	        return s;
	    } else if (v->data_type == MAT_T_UINT16) {
	        const uint16_t* p = static_cast<const uint16_t*>(v->data);
	        for (size_t i=0; i<N; ++i) 
	        	s.Cat((char)(p[i] & 0xFF)); // fallback simple (o hacer conversión real UTF-16→UTF-8)
	        return s;
	    }
	    // fallback UINT8
	    const uint8_t* p = static_cast<const uint8_t*>(v->data);
	    for (size_t i=0; i<N; ++i) 
	    	s.Cat((char)p[i]);
	    return s;
	}

    template <class T>
    static void ReadRowMajor(const matvar_t* v, Upp::Vector<T>& out, size_t& rows, size_t& cols) {
        Verify(v, "ReadMatrixRowMajor: null var");
        if constexpr (std::is_same<T,double>::value) 
        	Verify(v->class_type == MAT_C_DOUBLE, "type mismatch");
        if constexpr (std::is_same<T,int>::value) 
        	Verify(v->class_type == MAT_C_INT64, "type mismatch");
        Verify(v->rank == 2, "expected rank-2");
        rows = v->dims[0]; 
        cols = v->dims[1];
        const size_t n = rows*cols; 
        out.SetCount((int)n);
        Verify(v->data, "matrix has no data");
        
        auto row = ColToRow((const T*)v->data, rows, cols);
        std::memcpy(out.Begin(), row, n * sizeof(T));
    }
	template <class T>
	static Upp::Vector<T> ReadRowMajor(const matvar_t* v) {
		Upp::Vector<T> ret;
		size_t rows, cols;
		ReadRowMajor(v, ret, rows, cols);
		return ret;
	}
	    
	template <class T>
	static void ReadColMajor(const matvar_t* v, Upp::Vector<T>& out, size_t& rows, size_t& cols) {
	    Verify(v, "ReadMatrixColMajor: null var");
	    if constexpr (std::is_same<T, double>::value)
	        Verify(v->class_type == MAT_C_DOUBLE, "type mismatch");
	    if constexpr (std::is_same<T, int>::value)
	        Verify(v->class_type == MAT_C_INT64, "type mismatch");
	    Verify(v->rank == 2, "expected rank-2");
	    rows = v->dims[0];
	    cols = v->dims[1];
	    const size_t n = rows * cols;
	    out.SetCount((int)n);
	    Verify(v->data, "matrix has no data");
	
	    // Copy 1:1 in column-major order (same layout as MatIO buffer).
	    const T* src = static_cast<const T*>(v->data);
	    std::memcpy(out.Begin(), src, n * sizeof(T));
	}
	template <class T>
	static Upp::Vector<T> ReadColMajor(const matvar_t* v) {
		Upp::Vector<T> ret;
		size_t rows, cols;
		ReadColMajor(v, ret, rows, cols);
		return ret;
	}
		
    static void ReadRowMajor(const matvar_t* v, Upp::Vector<std::complex<double>>& out, size_t& rows, size_t& cols) {
        Verify(v && v->class_type == MAT_C_DOUBLE && v->isComplex, "ReadComplexMatrixRowMajor: not complex<double> matrix");
        Verify(v->rank == 2, "expected rank-2");
        rows = v->dims[0]; 
        cols = v->dims[1];
        const size_t n = rows*cols; 
        out.SetCount((int)n);
        const auto* cs = static_cast<const mat_complex_split_t*>(v->data);
        Verify(cs && cs->Re && cs->Im, "complex split missing");
        const double* re = (const double*)cs->Re;
        const double* im = (const double*)cs->Im;
        for (size_t c=0;c<cols;++c)
            for (size_t r=0;r<rows;++r)
                out[(int)(r*cols + c)] = std::complex<double>(re[c*rows + r], im[c*rows + r]);
    }
	static void ReadColMajor(const matvar_t* v, Upp::Vector<std::complex<double>>& out, size_t& rows, size_t& cols) {
	    Verify(v && v->class_type == MAT_C_DOUBLE && v->isComplex, "ReadComplexMatrixColMajor: not complex<double> matrix");
	    Verify(v->rank == 2, "expected rank-2");
	    rows = v->dims[0];
	    cols = v->dims[1];
	    const size_t n = rows * cols;
	    out.SetCount((int)n);
	    const auto* cs = static_cast<const mat_complex_split_t*>(v->data);
	    Verify(cs && cs->Re && cs->Im, "complex split missing");
	    const double* re = static_cast<const double*>(cs->Re);
	    const double* im = static_cast<const double*>(cs->Im);
	
	    for (size_t i = 0; i < n; ++i)		// MatIO stores re/im buffers in column-major order already.
	        out[(int)i] = std::complex<double>(re[i], im[i]);
	}


    static Eigen::MatrixXd ReadEigenDouble(const matvar_t* v) {
        Verify(v && v->class_type == MAT_C_DOUBLE && v->rank==2, "ReadEigenDouble: expected double matrix");
        Eigen::MatrixXd M(v->dims[0], v->dims[1]); // column-major
        Verify(v->data, "matrix has no data");
        std::memcpy(M.data(), v->data, (size_t)(v->dims[0]*v->dims[1]) * sizeof(double));
        return M;
    }
    static Eigen::MatrixXi ReadEigenInt(const matvar_t* v) {
        Verify(v && v->class_type == MAT_C_INT64 && v->rank==2, "ReadEigenInt: expected int matrix");
        Eigen::MatrixXi M(v->dims[0], v->dims[1]);
        Verify(v->data, "matrix has no data");
        std::memcpy(M.data(), v->data, (size_t)(v->dims[0]*v->dims[1]) * sizeof(int));
        return M;
    }
    static Eigen::MatrixXcd ReadEigenComplex(const matvar_t* v) {
        Verify(v && v->class_type == MAT_C_DOUBLE && v->isComplex && v->rank==2, "ReadEigenComplex: expected complex<double> matrix");
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> M(v->dims[0], v->dims[1]);
        const auto* cs = static_cast<const mat_complex_split_t*>(v->data);
        Verify(cs && cs->Re && cs->Im, "complex split missing");
        const double* re = (const double*)cs->Re; const double* im = (const double*)cs->Im;
        size_t rows = v->dims[0], cols = v->dims[1];
        for (size_t c=0;c<cols;++c)
            for (size_t r=0;r<rows;++r)
                M((Eigen::Index)r, (Eigen::Index)c) = std::complex<double>(re[c*rows + r], im[c*rows + r]);
        return M;
    }
    
	template <class T>
	static MultiDimMatrix<T> ReadMultiDim(const matvar_t* v) {
	    Verify(v, "ReadMultiDim: null matvar");
	    Verify(v->data, "ReadMultiDim: var has no data");
	
	    if constexpr (std::is_same<T, double>::value)
	        Verify(v->class_type == MAT_C_DOUBLE, "ReadMultiDim<double>: type mismatch");
	    else if constexpr (std::is_same<T, int>::value)
	        Verify(v->class_type == MAT_C_INT64, "ReadMultiDim<int>: type mismatch");
	    else if constexpr (std::is_same<T, std::complex<double>>::value)
	        Verify(v->class_type == MAT_C_DOUBLE && v->isComplex, "ReadMultiDim<complex>: expected complex<double> array");
	    else
	        static_assert(sizeof(T) == 0, "Unsupported MultiDimMatrix<T> (use double, int, std::complex<double>)");
	
	    Vector<int> dims;
	    dims.SetCount((int)v->rank);
	    for (int i = 0; i < dims.size(); ++i)
	        dims[i] = (int)v->dims[i];
	
	    MultiDimMatrix<T> M;
	    M.Resize(dims);
	    const size_t N = (size_t)M.size();
	
	    if constexpr (std::is_same<T, double>::value)
	        std::memcpy(M.begin(), v->data, N * sizeof(double));
	    else if constexpr (std::is_same<T, int>::value)
	        std::memcpy(M.begin(), v->data, N * sizeof(int));
	    else { // complex<double>
	        const auto* cs = static_cast<const mat_complex_split_t*>(v->data);
	        Verify(cs && cs->Re && cs->Im, "ReadMultiDim<complex>: split storage missing");
	        const double* re = (const double*)cs->Re;
	        const double* im = (const double*)cs->Im;
	        auto* out = M.begin();
	        for (size_t i = 0; i < N; ++i) out[i] = {re[i], im[i]};
	    }
	    return M;
	}

    void Close() noexcept { 
    	if (fp_) { 
    		Mat_Close(fp_); 
    		fp_ = nullptr; 
    	} 
    }
    
private:
    explicit MatFile(mat_t* fp, String p) : fp_(fp), path_(pick(p)) {}

    void Move(MatFile&& o) noexcept { 
    	fp_ = o.fp_; 
    	path_ = pick(o.path_); 
    	o.fp_ = nullptr; 
    }
	static inline Vector<size_t> ToSizeDims(const Vector<int>& vd) {
	    Vector<size_t> ds;
	    ds.SetCount(vd.size());
	    for (int i = 0; i < vd.size(); ++i)
	        ds[i] = (size_t)vd[i];
	    return ds;
	}
	
    mat_t* fp_ = nullptr;
    String path_;
};

} 

#endif