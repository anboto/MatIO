// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#ifndef _plugin_matio_matio_h_
#define _plugin_matio_matio_h_

#include <plugin/matio/lib/matio.h>
#include <Eigen/Eigen.h>


namespace Upp {

template<class T> void GetTypeCode 				(enum matio_classes &class_type, enum matio_types &data_type) {
	NEVER_("Class unsupported");
}

template<> void inline GetTypeCode<double>		(enum matio_classes &class_type, enum matio_types &data_type) {
	class_type = MAT_C_DOUBLE;		
	data_type  = MAT_T_DOUBLE;
}

template<> void inline GetTypeCode<float>		(enum matio_classes &class_type, enum matio_types &data_type) {
	class_type = MAT_C_SINGLE;		
	data_type  = MAT_T_SINGLE;
}

template<> void inline GetTypeCode<int>			(enum matio_classes &class_type, enum matio_types &data_type) {
	class_type = MAT_C_INT64;		
	data_type  = MAT_T_INT64;
}

template<> void inline GetTypeCode<String>		(enum matio_classes &class_type, enum matio_types &data_type) {
	class_type = MAT_C_CHAR;		
	data_type  = MAT_T_UTF8;
}

template<> void inline GetTypeCode<const char *> (enum matio_classes &class_type, enum matio_types &data_type) {
	class_type = MAT_C_CHAR;		
	data_type  = MAT_T_UTF8;
}


class MatVar {
public:
	MatVar() : var(nullptr), del(true) {}
	MatVar(const Nuller&)			{var = 0;} 
	bool IsNullInstance() const    	{return var == 0;}
	
	MatVar(String name, int rows, int cols, Vector<String> &elements) {
		Buffer<size_t> dim(2);
		dim[0] = rows;
		dim[1] = cols;
		
		Buffer<const char *> names(elements.GetCount());
		for (int i = 0; i < elements.GetCount(); ++i)
			names[i] = ~(elements[i]);
		var = Mat_VarCreateStruct(name, 2, dim, names, elements.GetCount());
		this->del = false;
	}
	
	~MatVar();
	
	const char *GetName()		{ASSERT(var != NULL); return var->name;}
	
	enum matio_classes GetType() 	{ASSERT(var != NULL); return var->class_type;}
	
	const char* GetTypeString();

	int GetDimCount() const			{ASSERT(var != NULL); return var->rank;}
	
	int GetDimCount(int dim) const	{ASSERT(var != NULL); return (int)var->dims[dim];}
	
	int GetCount();
	
	int GetFieldCount() {
		ASSERT(var != NULL); 
		ASSERT(var->class_type == MAT_C_STRUCT);
		
		return Mat_VarGetNumberOfFields(var);
	}
	const char *GetFieldName(int id) {	
		ASSERT(var != NULL); 
		ASSERT(var->class_type == MAT_C_STRUCT);
		
		return Mat_VarGetStructFieldnames(var)[id];
	}
	
	MatVar GetVar(String name) {
		ASSERT(var != NULL); 
		ASSERT(var->class_type == MAT_C_STRUCT);
		
		return MatVar(Mat_VarGetStructFieldByName(var, name, 0), false);
	}
	
	bool IsLoaded() {return var != 0;}
	
	
	void SetStruct(String name, MatVar &val) {
		Mat_VarSetStructFieldByName(var, name, 0, val.var); 
	}
	
	template <class T>
	void SetStruct(String name, const void *data, int numRows, int numCols, int index) {
		int numDim = 2;
		if (IsNull(numCols)) 
			numDim = 1;
			
		Buffer<size_t> dims(numDim);
		dims[0] = numRows;
		if (!IsNull(numCols)) 
			dims[1] = numCols;
		
		enum matio_classes class_type;
		enum matio_types data_type;
		GetTypeCode<T>(class_type, data_type);
		
		matvar_t *variable = Mat_VarCreate(name, class_type, data_type, numDim, dims, data, MAT_F_DONT_COPY_DATA);
		Mat_VarSetStructFieldByName(var, name, index, variable); 
	}

	void SetStruct(String name, Eigen::MatrixXd &data, int index = 0) {
		SetStruct<double>(name, data.data(), (int)data.rows(), (int)data.cols(), index);
	}
	
	void SetStruct(String name, String data, int index = 0) {
		SetStruct<String>(name, (void *)data.Begin(), 1, data.GetCount(), index);
	}

	void SetStruct(String name, const char *data, int index = 0) {
		SetStruct<String>(name, (void *)data, 1, (int)strlen(data), index);
	}

	template <class T>
	void SetStruct(String name, T data, int index = 0) {
		SetStruct<T>(name, &data, 1, Null, index);
	}
	
private:
	MatVar(mat_t *mat, String name);
	MatVar(matvar_t *var, bool del) {this->var = var; this->del = del;}
	
	matvar_t *var = nullptr;
	bool del;
	
	friend class MatFile;
};
	
class MatFile {
public:
	MatFile() : mat(NULL), listVar(NULL), numVar(0) {Mat_LogInitFunc("Matio", LogFunc);}
	~MatFile();
	
	bool Create(String fileName, mat_ft version = MAT_FT_MAT5);
	void Close();
	
	bool OpenRead(String fileName) 	{return Open(fileName, MAT_ACC_RDONLY);}
	bool OpenWrite(String fileName) {return Open(fileName, MAT_ACC_RDWR);}
	bool IsOpen()					{return !!mat;}
	bool OpenCreate(String fileName, enum mat_ft mat_file_ver) {
		if (mat != NULL)
			Mat_Close(mat);
	
		mat = Mat_CreateVer(fileName, "BEMRosetta", mat_file_ver);
			    
    	return !!mat;
	}
	
	mat_ft GetVersion() {
		ASSERT(mat != NULL);
		
		return Mat_GetVersion(mat);
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
	
	int GetVarCount() {
		GetVarList();
		return (int)numVar;
	}
	
	String GetVarName(int id) {
		GetVarList();
		if (id >= (int)numVar)
			return Null;
		return listVar[id];
	}
		
	bool Exist(String name, bool nocase = false) {
		GetVarList();
		if (numVar == 0)
			return false;
		
		if (nocase) {
			name = ToLower(name);
			name.Replace("_", "");
		}
		
		for (int i = 0; i < (int)numVar; ++i) {
			if (nocase) {
				String var = ToLower(listVar[i]);
				var.Replace("_", "");
				if (var == name)
					return true;
			} else {
				if (String(listVar[i]) == name)
					return true;
			}
		}
		return false;
	}
	
	int Exist(const Vector<String> &name, bool nocase = false) {
		GetVarList();
		if (numVar == 0)
			return -1;
		
		Vector<String> names;
		if (!nocase)
			names = clone(name);
		else {
			names.SetCount(name.size());
			for (int i = 0; i < name.size(); ++i) {
				names[i] = ToLower(name[i]);
				names[i].Replace("_", "");
			}
		}
		
		for (int i = 0; i < (int)numVar; ++i) {
			String var;
			if (nocase) {
				var = ToLower(listVar[i]);
				var.Replace("_", "");
			} else
				var = String(listVar[i]);
			for (int iv = 0; iv < name.size(); ++iv) {
				if (var == names[iv])
					return iv;
			}
		}
		return -1;
	}
	
	bool VarDelete(String name) {
		ASSERT(mat != NULL);
		
		int ret = Mat_VarDelete(mat, name);
		GetVarList();
		if (ret == 0) 
			return true;
		return false;
	}
	
	MatVar GetVar(String name) const {return MatVar(mat, name);}

	template <class T> inline
	T Get(const MatVar &var) {
		ASSERT(mat != NULL);
	
		T ret = Null;
			
		int numDim = var.GetDimCount();	
		if (numDim > 2)
			return ret;
			
		Buffer<int> start(numDim, 0);
		Buffer<int> stride(numDim, 1);
		Buffer<int> edge(numDim);
	
		for (int i = 0; i < numDim; ++i)
			edge[i] = var.GetDimCount(i);
	
		if (0 != Mat_VarReadData(mat, var.var, &ret, start, stride, edge)) 
			return ret;
		
		return ret;	
	}
	
	MatVar GetVarThrow(String name, bool nocase = false) {
		if (nocase) {
			if (numVar == 0)
				throw Exc(t_("Matio: No vars"));
		
			name = ToLower(name);
			name.Replace("_", "");
			int i;
			for (i = 0; i < (int)numVar; ++i) {
				String var = ToLower(listVar[i]);
				var.Replace("_", "");
				if (var == name) {
					name = listVar[i];
					break;
				}
			}
			if (i == (int)numVar)
				throw Exc(t_("Matio: Var does not exist"));
		}
		return GetVar(name);
	}
	
	template <class T> inline
	T Get(String name, bool nocase = false) {
		MatVar var = GetVarThrow(name, nocase);	
		return Get<T>(var);
	}
	
	template <class T>
	void Get(MatVar &var, Vector<T> &ret) {
		ASSERT(mat != NULL);
			
		int numDim = var.GetDimCount();	
		if (numDim > 2)
			throw Exc("Matio: Size is not 2");
			
		Buffer<int> start(numDim, 0);
		Buffer<int> stride(numDim, 1);
		Buffer<int> edge(numDim);
	
		int sz = 1;
		for (int i = 0; i < numDim; ++i) {
			edge[i] = var.GetDimCount(i);
			sz *= var.GetDimCount(i);
		}
		ret.SetCount(sz);
	
		if (0 != Mat_VarReadData(mat, var.var, ret.begin(), start, stride, edge)) 
			throw Exc("Matio: Problem reading var");
	}

	template <class T>
	void Get(String name, Vector<T> &ret, bool nocase = false) {
		MatVar var = GetVarThrow(name, nocase);
		Get<T>(var, ret);
	}
	
	void Get(MatVar &var, Eigen::VectorXd &ret) {
		ASSERT(mat != NULL);
			
		int numDim = var.GetDimCount();	
		if (numDim > 2)
			throw Exc("Matio: Size is not 2");
			
		Buffer<int> start(numDim, 0);
		Buffer<int> stride(numDim, 1);
		Buffer<int> edge(numDim);
	
		for (int i = 0; i < numDim; ++i)
			edge[i] = var.GetDimCount(i);
	
		ret.resize((int)(var.GetDimCount(0)*(int)var.GetDimCount(1)));
	
		if (0 != Mat_VarReadData(mat, var.var, ret.data(), start, stride, edge)) 
			throw Exc("Matio: Problem reading var");
	}

	void Get(String name, Eigen::VectorXd &ret, bool nocase = false) {
		MatVar var = GetVarThrow(name, nocase);
		Get(var, ret);
	}
	
	void Get(MatVar &var, Eigen::MatrixXd &ret) {
		ASSERT(mat != NULL);
			
		int numDim = var.GetDimCount();	
		if (numDim > 2)
			throw Exc("Matio: Size is not 2");
			
		Buffer<int> start(numDim, 0);
		Buffer<int> stride(numDim, 1);
		Buffer<int> edge(numDim);
	
		for (int i = 0; i < numDim; ++i)
			edge[i] = var.GetDimCount(i);
	
		//Buffer<double> d((int)var.GetDimCount(0)*(int)var.GetDimCount(1));
		ret.resize(var.GetDimCount(0), var.GetDimCount(1));
		if (0 != Mat_VarReadData(mat, var.var, ret.data(), start, stride, edge)) 
			throw Exc("Matio: Problem reading var");
		
		//CopyRowMajor(d.Get(), int(var.GetDimCount(0)), int(var.GetDimCount(1)), ret);
	}
	
	void Get(String name, Eigen::MatrixXd &ret, bool nocase = false) {
		MatVar var = GetVarThrow(name, nocase);
		Get(var, ret);
	}
	
	template <class T>
	void Get(MatVar &var, MultiDimMatrix<std::complex<T>> &ret) {
		ASSERT(mat != NULL);
			
		int numDim = var.GetDimCount();	
			
		Buffer<int> start(numDim, 0);
		Buffer<int> stride(numDim, 1);
		Buffer<int> edge(numDim);
		Vector<int> vedge(numDim);
		
		for (int i = 0; i < numDim; ++i)
			edge[i] = vedge[i] = var.GetDimCount(i);
	
		Buffer<double> real(var.GetCount()), imag(var.GetCount());
	
		struct mat_complex_split_t data = {real, imag};
		
		if (0 != Mat_VarReadData(mat, var.var, &data, start, stride, edge)) 
			throw Exc("Matio: Problem reading var");
		
		ret.Resize(vedge);
		for (int i = 0; i < ret.size(); ++i) {
			ret.begin()[i].real(real[i]);
			ret.begin()[i].imag(imag[i]);
		}		
	}	

	template <class T>
	void Get(String name, MultiDimMatrix<std::complex<T>> &ret, bool nocase = false) {
		MatVar var = GetVarThrow(name, nocase);
		Get<T>(var, ret);
	}
	
	template <class T>
	typename std::enable_if<std::is_floating_point_v<T>>::type
	Get(MatVar &var, MultiDimMatrix<T> &ret) {
		ASSERT(mat != NULL);
			
		int numDim = var.GetDimCount();	
			
		Buffer<int> start(numDim, 0);
		Buffer<int> stride(numDim, 1);
		Buffer<int> edge(numDim);
		Vector<int> vedge(numDim);
		
		for (int i = 0; i < numDim; ++i)
			edge[i] = vedge[i] = var.GetDimCount(i);
	
		ret.Resize(vedge);
	
		if (0 != Mat_VarReadData(mat, var.var, ret.begin(), start, stride, edge))
			throw Exc("Matio: Problem reading var");
	}	

	template <class T>
	void Get(String name, MultiDimMatrix<T> &ret, bool nocase = false) {
		MatVar var = GetVarThrow(name, nocase);
		Get<T>(var, ret);
	}
			
	void Get(MatVar &var, Eigen::MatrixXcd &ret) {
		ASSERT(mat != NULL);
			
		int numDim = var.GetDimCount();	
		if (numDim > 2)
			throw Exc("Matio: Size is not 2");
			
		Buffer<int> start(numDim, 0);
		Buffer<int> stride(numDim, 1);
		Buffer<int> edge(numDim);
	
		for (int i = 0; i < numDim; ++i)
			edge[i] = var.GetDimCount(i);
	
		Buffer<double> real(var.GetCount()), imag(var.GetCount());
		
		struct mat_complex_split_t data = {real, imag};
		
		if (0 != Mat_VarReadData(mat, var.var, &data, start, stride, edge)) 
			throw Exc("Matio: Problem reading var");

		ret.resize(edge[0], edge[1]);		
		for (int i = 0; i < ret.size(); ++i) {
			ret(i).real(real[i]);
			ret(i).imag(imag[i]);
		}
	}
	
	void Get(String name, Eigen::MatrixXcd &ret, bool nocase = false) {
		MatVar var = GetVarThrow(name, nocase);
		Get(var, ret);
	}
	
	void GetCell(String name, Vector<String> &ret) {
		ret.Clear();
	    matvar_t *cell_array = Mat_VarRead(mat, name);
	    if (!cell_array) 
	        throw Exc(Format(t_("Variable '%s' not found"), name));
	      
	    if (cell_array->class_type != MAT_C_CELL) {
	        Mat_VarFree(cell_array);
	        throw Exc(t_("Variable is not a cell array"));
	    }
	    size_t num = 1;
	    for (int i = 0; i < cell_array->rank; i++) 
	        num *= cell_array->dims[i];
	    
	    for (size_t i = 0; i < num; i++) {
	        matvar_t *element = ((matvar_t **)cell_array->data)[i];
	        if (!element)
	            continue;
	
	        if (element->class_type == MAT_C_CHAR) {
	            int total_chars = int(element->dims[0] * element->dims[1]);
	
	            StringBuffer c_str(total_chars + 1);
				if (element->data_type == MAT_T_UINT16) {
		            uint16_t *matlab_str = (uint16_t*)element->data;
		            for (int j = 0; j < total_chars; j++)
		                c_str[j] = (char)matlab_str[j];
				} else if (element->data_type == MAT_T_UTF8)
		            memcpy(c_str, element->data, total_chars);
				
				c_str[total_chars] = '\0';
				ret << String(c_str.begin(), total_chars);
	        } else 
	        	throw Exc(Format(t_("Cell element %d is not a string"), (int)i));
	    }
	    Mat_VarFree(cell_array);
	}
	
	void SetCell(String name, const Vector<String> &strs, bool compression = false) {
		size_t num = strs.size();
		Buffer<matvar_t *> cell_elements(num);
	
	    for (int i = 0; i < (int)num; i++) {
	        const char *c_str = strs[i];
	
	        size_t dims[2] = {1, strlen(c_str)};
	
	        cell_elements[i] = Mat_VarCreate(NULL, MAT_C_CHAR, MAT_T_UTF8, 2, dims, c_str, MAT_F_DONT_COPY_DATA);
	    }
	    size_t cell_dims[2] = {num, 1};
	
	    matvar_t *cell_array = Mat_VarCreate(name, MAT_C_CELL, MAT_T_CELL, 2, cell_dims, cell_elements, MAT_F_DONT_COPY_DATA);
	
	    Mat_VarWrite(mat, cell_array, compression ? MAT_COMPRESSION_NONE : MAT_COMPRESSION_ZLIB);
	
	    for (size_t i = 0; i < num; i++) 
	        Mat_VarFree(cell_elements[i]);
	    
	    Mat_VarFree(cell_array);
	}	
	
	bool Set(MatVar &var, bool compression = true) {
		if (0 != Mat_VarWrite(mat, var.var, compression ? MAT_COMPRESSION_NONE : MAT_COMPRESSION_ZLIB))
			return false;
		return true;
	}

	template <class T>
	bool Set(String name, const void *data, const Buffer<size_t> &dims, int numDim, bool compression = true) {
		enum matio_classes class_type;
		enum matio_types data_type;
		GetTypeCode<T>(class_type, data_type);
		
		if (Exist(name))
			if (0 != Mat_VarDelete(mat, name))
				return false;
		matvar_t *var = Mat_VarCreate(name, class_type, data_type, numDim, dims, data, MAT_F_DONT_COPY_DATA);
		if (var == NULL)
			return false;
		if (0 != Mat_VarWrite(mat, var, compression ? MAT_COMPRESSION_NONE : MAT_COMPRESSION_ZLIB))
			return false;
		return true;
	}
	
	template <class T>
	bool Set(String name, const void *data, const Vector<int> &vdims, bool compression = true) {
		Buffer<size_t> dims(vdims.size());
		for (int i = 0; i < vdims.size(); ++i)
			dims[i] = vdims[i];
		
		return Set<T>(name, data, dims, vdims.size(), compression);
	}
			
	template <class T>
	bool Set(String name, const void *data, int numRows, int numCols, bool compression = true) {
		int numDim = 2;
		if (IsNull(numCols)) 
			numDim = 1;
			
		Buffer<size_t> dims(numDim);
		dims[0] = numRows;
		if (!IsNull(numCols)) 
			dims[1] = numCols;
		
		return Set<T>(name, data, dims, numDim, compression);
	}
	
	template <class T>
	bool Set(String name, const Vector<T> &data, bool compression = true) {
		return Set<T>(name, data, data.size(), Null, compression);
	}

	bool Set(String name, const Eigen::VectorXd &data, bool compression = true) {
		return Set<double>(name, data.data(), (int)data.size(), Null, compression);
	}
	
	bool Set(String name, const Eigen::MatrixXd &data, bool compression = true) {
		return Set<double>(name, data.data(), (int)data.rows(), (int)data.cols(), compression);
	}

	template <class T>
	bool Set(String name, const MultiDimMatrix<std::complex<T>> &data, bool compression = true) {
		Buffer<size_t> dims(data.GetNumAxis());
		for (int i = 0; i < data.GetNumAxis(); ++i)
			dims[i] = data.GetAxisDim(i);
		return SetComplex<T>(name, data.begin(), dims, data.GetNumAxis(), compression);
	}
	
	template <class T>
	bool Set(String name, const MultiDimMatrix<T> &data, bool compression = true) {
		Buffer<size_t> dims(data.GetNumAxis());
		for (int i = 0; i < data.GetNumAxis(); ++i)
			dims[i] = data.GetAxisDim(i);
		return Set<T>(name, data.begin(), dims, data.GetNumAxis(), compression);
	}
	
	bool Set(String name, String data, bool compression = true) {
		return Set<String>(name, (void *)data.Begin(), 1, data.GetCount(), compression);
	}
	
	bool Set(String name, const char *data, bool compression = true) {
		return Set<String>(name, data, 1, (int)strlen(data), compression);
	}
			
	template <class T>
	bool Set(String name, T data, bool compression = true) {
		return Set<T>(name, &data, 1, Null, compression);
	}
	
	template <class T>
	bool SetComplex(String name, const void *data, const Buffer<size_t> &dims, int numDim, bool compression = true) {
		enum matio_classes class_type;
		enum matio_types data_type;
		GetTypeCode<T>(class_type, data_type);
		
		if (Exist(name))
			if (0 != Mat_VarDelete(mat, name))
				return false;
			
		int sz = 1;
		for (int i = 0; i < numDim; ++i)
			sz *= (int)dims[i];
		Buffer<double> real(sz), imag(sz);
		std::complex<T> *cdata = (std::complex<T>*)data;
		for (int i = 0; i < sz; ++i) {
			real[i] = cdata[i].real();
			imag[i] = cdata[i].imag();
		}
		
		struct mat_complex_split_t sdata = {real, imag};
		
		matvar_t *var = Mat_VarCreate(name, class_type, data_type, numDim, dims, &sdata, MAT_F_DONT_COPY_DATA | MAT_F_COMPLEX);
		if (var == NULL)
			return false;
		if (0 != Mat_VarWrite(mat, var, compression ? MAT_COMPRESSION_NONE : MAT_COMPRESSION_ZLIB))
			return false;
		return true;
	}
	
	bool Set(String name, const Eigen::MatrixXcd &data, bool compression = true) {
		Buffer<size_t> dims(2);
		dims[0] = data.rows();
		dims[1] = data.cols();
		return SetComplex<double>(name, data.data(), dims, 2, compression);
	}
		
	mat_t *mat = nullptr;
	
private:
	char *const* listVar = nullptr;
	size_t numVar;
	
	bool Open(String fileName, int mode) {
		if (mat != NULL)
			Mat_Close(mat);
	
		mat = Mat_Open(fileName, mode);
			    
    	return !!mat;
	}
	
	void GetVarList() {
		ASSERT(mat != NULL);
		
		numVar = 0; 
		listVar = Mat_GetDir(mat, &numVar);
	}
	
	static void LogFunc(int log_level, const char *message) {
		throw Exc(message);
	}
};

template <> inline
String MatFile::Get<String>(const MatVar &var) {
	ASSERT(mat != NULL);
		
	int numDim = var.GetDimCount();	
	if (numDim > 2)
		return String();
		
	Buffer<int> start(numDim, 0);
	Buffer<int> stride(numDim, 1);
	Buffer<int> edge(numDim);

	for (int i = 0; i < numDim; ++i)
		edge[i] = var.GetDimCount(i);

	StringBuffer ret(var.GetDimCount(0)*var.GetDimCount(1));

	var.var->class_type = MAT_C_INT8;
	var.var->data_type  = MAT_T_INT8;

	if (0 != Mat_VarReadData(mat, var.var, (void *)ret.Begin(), start, stride, edge)) 
		return String();
	
	return String(ret);	
}

}
	
#endif

