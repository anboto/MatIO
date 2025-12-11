// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2025, the Anboto author and contributors
#include <Core/Core.h>
#include <Functions4U/Functions4U.h>
#include <MatIO/matio.h>
#include "data.brc"


using namespace Upp;

static void ReadDemo() {
	String fname = "data.mat";
	
	VERIFY(SaveFile(fname, String(dataMat, dataMat_length))); 
	
	MatFile f;
	VERIFY(f.OpenRead(fname));

	VERIFY(f.Exist("H"));
	MatVar hVar = f.GetVar("H");
    VERIFY(abs(MatFile::ReadDouble(hVar) - 10) < 1e-12);
    
	MatVar Bvar = f.GetVar("B");
	Eigen::VectorXd B = MatFile::ReadEigenDouble(Bvar);
	VERIFY(B.rows()==100 && B.cols()==1);
	VERIFY(abs(B(0)-8.411362) < 0.001 && abs(B(1)-66.79221) < 0.001 && abs(B(2)-222.6293) < 0.001);
	
	MatVar dataVar = f.GetVar("Data");
	MatVar nameVar = dataVar.GetStructVar("name");
	VERIFY(nameVar);
	String str = MatFile::ReadString(nameVar);
	VERIFY(str == "Speed");
	
	MatVar optionsVar = f.GetVar("Options");
	MatVar frecVar = optionsVar.GetStructVar("FreqChoice");
	VERIFY(frecVar);
	String strf = MatFile::ReadString(frecVar);
	VERIFY(strf == "G");
	
	MatVar optimVar = optionsVar.GetStructVar("Optim");
	VERIFY(optimVar);
	VERIFY(optimVar.ExistStructVar("maxEval"));
	MatVar maxVar = optimVar.GetStructVar("maxEval");
	VERIFY(maxVar);
	VERIFY(MatFile::ReadInt(maxVar) == 100);
	
	Cout() << "ReadDemo OK\n";
}

static void DemoStruct(mat_ft ver, const String& fname) {
    {
        MatFile f;
        f.OpenCreate(fname, ver);

        f.Write("answer", (int)42);
        f.Write("pi", 3.141592653589793);
        f.Write("z1", std::complex<double>(2.5, -1.25));
        f.Write("greeting", "Kaixo Iñaki!");

        {	// Matrix (2x3) column-major data: [1 3 5; 2 4 6]
            const size_t R = 2, C = 3;
            double Acol[] = {1,2, 3,4, 5,6};
            f.WriteColMajor("A", R, C, Acol);
        }
        {	// Int matrix (2x3) row-major data
            const size_t R = 2, C = 3;
            int Mrow[] = {10,11,12, 20,21,22};
            f.WriteRowMajor("M", R, C, Mrow);
        }
        {	// Complex (2x2) matrix, provided row-major then converted internally
            std::complex<double> data[] = {std::complex<double>(1,0.5),  std::complex<double>(2,-0.5),
                  						   std::complex<double>(3,1.0),  std::complex<double>(4,-1.0)};
            f.WriteRowMajor("Cplx", 2, 2, data);
        }
		{	// Eigen (double, int, complex<double>)
            Eigen::MatrixXd E(2,3); E << 1,3,5, 2,4,6;         // col-major by default
            f.Write("E", E);

            Eigen::MatrixXi Ei(2,3);
            Ei << 10,11,12, 20,21,22;
            f.Write("Ei", Ei);

            Eigen::MatrixXcd Ec(2,2);
            Ec << std::complex<double>(1,0.5),  std::complex<double>(2,-0.5),
                  std::complex<double>(3,1.0),  std::complex<double>(4,-1.0);
            f.Write("Ec", Ec);
        }
		{
			Upp::MultiDimMatrix<double> A_md(4, 3, 2); // dims {4,3,2}, column-major
			A_md.SetConstant(1.23);
			A_md(1, 2, 1) = 3.21;
			f.Write("A_md", A_md);
		}

        MatFile::StructNode options(f, "Options", {"Optim", "project", "Z"});

        MatFile::StructNode& optim = options.AddChild("Optim", {"maxEval","tol","method","zstart","H"});
        optim.Write("maxEval", (int)1000);
        optim.Write("tol", 1e-8);
        optim.Write("method", "BFGS");
        optim.Write("zstart", std::complex<double>(0.0, 1.0));
        {	// 2x2 complex H in column-major (matches how MATLAB stores)
            std::complex<double> Cm[] = {std::complex<double>(1,0.5),  std::complex<double>(2,-0.5),
                  						 std::complex<double>(3,1.0),  std::complex<double>(4,-1.0)};
            optim.WriteColMajor("H", 2, 2, Cm);
        }
        options.Write("project", "Hydrodynamic calculation");
		{
			Upp::MultiDimMatrix<std::complex<double>> Z(5, 4, 3);
			Z.SetZero(); Z(1,2,1) = {2.0, -1.0};
			options.Write("Z", Z);  // overload above
		}
        options.Write(MAT_COMPRESSION_ZLIB);
    }
    {
        MatFile f;
        f.OpenRead(fname);

        // Scalars
        {
            MatVar piVar = f.GetVar("pi");
            VERIFY(std::abs(MatFile::ReadDouble(piVar) - 3.141592653589793) < 1e-12);
        }
        {
            MatVar ansVar = f.GetVar("answer");
            VERIFY(MatFile::ReadInt(ansVar) == 42);
        }
        {
            MatVar z1Var = f.GetVar("z1");
            std::complex<double> z = MatFile::ReadComplex(z1Var);
            VERIFY(std::abs(z.real() - 2.5) < 1e-12 && std::abs(z.imag() + 1.25) < 1e-12);
        }
        {
            MatVar grVar = f.GetVar("greeting");
            String s = MatFile::ReadString(grVar);
            VERIFY(s == "Kaixo Iñaki!");
        }
        
        {	// Matrices (double/int)
            MatVar Avar = f.GetVar("A");
            size_t r=0,c=0; 
            Vector<double> Arow;
            MatFile::ReadRowMajor<double>(Avar, Arow, r, c);
            VERIFY(r==2 && c==3);
            // Expect row-major [1 3 5 ; 2 4 6]
            VERIFY(Arow[0]==1 && Arow[1]==3 && Arow[2]==5 && Arow[3]==2 && Arow[4]==4 && Arow[5]==6);
        }
        {
            MatVar Mvar = f.GetVar("M");
            size_t r=0,c=0; 
            Vector<int> Mrow2;
            MatFile::ReadRowMajor<int>(Mvar, Mrow2, r, c);
            VERIFY(r==2 && c==3);
            VERIFY(Mrow2[0]==10 && Mrow2[1]==11 && Mrow2[2]==12 && Mrow2[3]==20 && Mrow2[4]==21 && Mrow2[5]==22);
        }
        {	// Complex matrix
            MatVar Cvar = f.GetVar("Cplx");
            size_t r=0,c=0; 
            Vector<std::complex<double>> Crow;
            MatFile::ReadRowMajor(Cvar, Crow, r, c);
            VERIFY(r==2 && c==2);
            VERIFY(abs(Crow[0].real()-1.0) < 1e-12 && abs(Crow[0].imag()-0.5) < 1e-12);
        }

        {	// Eigen reads
            MatVar Evar = f.GetVar("E");
            Eigen::MatrixXd E = MatFile::ReadEigenDouble(Evar);
            VERIFY(E.rows()==2 && E.cols()==3);
            VERIFY(E(0,0)==1 && E(0,1)==3 && E(0,2)==5);
        }
        {
            MatVar Eivar = f.GetVar("Ei");
            Eigen::MatrixXi Ei = MatFile::ReadEigenInt(Eivar);
            VERIFY(Ei(0,0)==10 && Ei(1,2)==22);
        }
        {
            MatVar Ecvar = f.GetVar("Ec");
            Eigen::MatrixXcd Ec = MatFile::ReadEigenComplex(Ecvar);
            VERIFY(abs(Ec(0,0).real()-1.0) < 1e-12 && abs(Ec(0,0).imag()-0.5) < 1e-12);
        }
		{			
			Upp::MultiDimMatrix<double> A_md = Upp::MatFile::ReadMultiDim<double>(f.GetVar("A_md"));
			VERIFY(A_md.size() == 24);        // 4*3*2
			VERIFY(A_md.GetAxisDim()[1] == 3);  // {4,3,2}
			VERIFY(abs(A_md(1, 1, 1) - 1.23) < 1E-12);
			VERIFY(abs(A_md(1, 2, 1) - 3.21) < 1E-12);
		}
        {	// Nested struct: Options → Optim → fields
            MatVar OptVar = f.GetVar("Options");
            MatVar Optim = OptVar.GetStructVar("Optim");
            VERIFY(Optim);

            MatVar maxEval = Optim.GetStructVar("maxEval");
            MatVar tol     = Optim.GetStructVar("tol");
            MatVar zstart  = Optim.GetStructVar("zstart");
            MatVar H       = Optim.GetStructVar("H");

            VERIFY(MatFile::ReadInt(maxEval) == 1000);
            VERIFY(abs(MatFile::ReadDouble(tol) - 1e-8) < 1e-15);

            std::complex<double> z0 = MatFile::ReadComplex(zstart);
            VERIFY(abs(z0.real()) < 1e-12 && abs(z0.imag()-1.0) < 1e-12);

            size_t Hr=0,Hc=0; 
            Vector<std::complex<double>> Hrow;
            MatFile::ReadRowMajor(H, Hrow, Hr, Hc);
            VERIFY(Hr==2 && Hc==2);
            
            VERIFY(MatFile::ReadString(OptVar.GetStructVar("project")) == "Hydrodynamic calculation");
            
            Upp::MultiDimMatrix<std::complex<double>> Z = Upp::MatFile::ReadMultiDim<std::complex<double>>(OptVar.GetStructVar("Z"));
            VERIFY(Z(0,0,0).real() == 0 && abs(Z(1,2,1).real() - 2) < 1E-12 && abs(Z(1,2,1).imag() + 1) < 1E-12);
        }
    }

    Cout() << "DemoStruct OK for " << fname << "\n";
}


CONSOLE_APP_MAIN
{
	Cout() << "\nMatIO demo\n";
	
	try {
		ReadDemo();
        DemoStruct(mat_ft::MAT_FT_MAT5,  "demo_v5.mat");
        DemoStruct(mat_ft::MAT_FT_MAT73, "demov73.mat");
	} catch(Exc e) {
		Cout() << "\nError: " << e;
		SetExitCode(1);
	}

	#ifdef flagDEBUG
	Cout() << "\nPress enter to end";
	ReadStdIn();
	#endif
}
