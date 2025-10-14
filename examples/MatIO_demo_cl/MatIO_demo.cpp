// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#include <Core/Core.h>
#include "data.brc"
#include <Functions4U/Functions4U.h>

using namespace Upp;

#include <MatIO/matio.h>


void PrintVariables(MatFile &mfile) {
	Cout() << "\nVariables list:";
	for (int iv = 0; iv < mfile.GetVarCount(); ++iv) {
		String name = mfile.GetVarName(iv);
		Cout() << Format("\n- %s", name);
		if (name != "#subsystem#" && name != "#refs#") {
			MatVar var = mfile.GetVar(name);
			Cout() << Format(" <%s>", var.GetTypeString());
			for (int id = 0; id < var.GetDimCount(); ++id) {
				if (id == 0)
					Cout() << "(";
				else
					Cout() << ",";
				Cout() << FormatInt(var.GetDimCount(id));
				if (id == var.GetDimCount() -1)
					Cout() << ")";
			}
		}
	}	
}

void DemoStructure() {
	MatFile mfile;
	
	String fileName = "demoStruct.mat";
	
	if (!FileExists(fileName)) 
		VERIFY(SaveFile(fileName, String(dataMat, dataMat_length))); // Error: Impossible to save .mat file
	
	VERIFY(mfile.OpenWrite(fileName));
		
	Cout() << "\nVersion " << mfile.GetVersionName();
	PrintVariables(mfile) ;
	
	MatVar str = mfile.GetVar("Structure");
	if (!str.IsLoaded()) 
		Cout() << "\nStructure was deleted";
	else {
		Cout() << "\nStructure " << str.GetTypeString();
		for (int i = 0; i < str.GetFieldCount(); ++i) {
			String name = str.GetFieldName(i);
			MatVar var = str.GetVar(name);
			Cout() << Format("\n- %s <%s>(%d, %d)", name, var.GetTypeString(), var.GetDimCount(0), var.GetDimCount(1));
		}
		//Cout() << "\nStructure.Est2.One: " << mfile.Get<int>(mfile.GetVar("Structure").GetVar("Est2").GetVar("One"));
		//Cout() << "\nStructure.Est2.Two: " << mfile.Get<String>(mfile.GetVar("Structure").GetVar("Est2").GetVar("Two"));
		
		mfile.VarDelete("Structure");
	}
	
	Vector<String> elements;
	elements << "First" << "Est2" << "Third";
	MatVar structure("Structure", 1, 1, elements);
	structure.SetStruct("First", "The first");

	Vector<String> elements2;
	elements2 << "One" << "Two";
	MatVar est2("Est2", 1, 1, elements2);
	est2.SetStruct("Two", "The second 22");
	est2.SetStruct("One", 211111);
	
	structure.SetStruct("Est2", est2);
	
	mfile.Set(structure);	
	
	MatVar opt = mfile.GetVar("Options");
	if (!opt.IsLoaded()) 
		Cout() << "\nOptions was deleted";
	else {
		Cout() << "\nOptions " << opt.GetTypeString();
		int numf = opt.GetFieldCount();
		for (int i = 0; i < numf; ++i) {
			String name = opt.GetFieldName(i);
			MatVar var = opt.GetVar(name);
			Cout() << Format("\n- %s <%s>(%d, %d)", name, var.GetTypeString(), var.GetDimCount(0), var.GetDimCount(1));
		}
		
		Cout() << "\nFreqRangeChoice: " << mfile.Get<String>(mfile.GetVar("Options").GetVar("FreqRangeChoice"));
		Cout() << "\nFreqChoice: "		<< mfile.Get<String>(mfile.GetVar("Options").GetVar("FreqChoice"));
		Cout() << "\nMode: " 			<< mfile.Get<double> (mfile.GetVar("Options").GetVar("Mode"));
		Cout() << "\nMod: " 			<< mfile.Get<double> (mfile.GetVar("Options").GetVar("Mod"));
	
	
		MatVar optim = opt.GetVar("Optim");
		int numOpt = optim.GetFieldCount();
		for (int i = 0; i < numOpt; ++i) {
			String name = optim.GetFieldName(i);
			MatVar var = optim.GetVar(name);
			Cout() << Format("\n- %s <%s>(%d, %d)", name, var.GetTypeString(), var.GetDimCount(0), var.GetDimCount(1));
		}
		
		Cout() << "\nInit.Cond: " 	  << mfile.Get<double>(mfile.GetVar("Options").GetVar("Optim").GetVar("InitCond"));
		Cout() << "\nInit.maxEval: "  << mfile.Get<double>(mfile.GetVar("Options").GetVar("Optim").GetVar("maxEval"));
		
		Cout() << "\nInit.Tol: " 	  << mfile.Get<double>(mfile.GetVar("Options").GetVar("Optim").GetVar("Tol"));
		Cout() << "\nInit.StepTol: "  << mfile.Get<double>(mfile.GetVar("Options").GetVar("Optim").GetVar("StepTol"));
		Cout() << "\nInit.ThresRel: " << mfile.Get<double>(mfile.GetVar("Options").GetVar("Optim").GetVar("ThresRel"));
		Cout() << "\nInit.ThresAbs: " << mfile.Get<double>(mfile.GetVar("Options").GetVar("Optim").GetVar("ThresAbs"));
	
		mfile.VarDelete("Options");
	}

	//Create a 1 x 1 struct 'Data' with fields: name, unit, value
	const char *structname = "Data";
	if (!mfile.Exist(structname)) {
		const char *fieldnames[3] = { "name","unit","value" };
		size_t structdim[2] = { 1, 1 }; // create 1 x 1 struct
		matvar_t* matstruct = Mat_VarCreateStruct(structname, 2, structdim, fieldnames, 3); 
	
		const char* mystring = "Speed";
		size_t dim[2] = { 1, 5 }; //string dimension
		matvar_t *variable = Mat_VarCreate(fieldnames[0], MAT_C_CHAR, MAT_T_UTF8, 2, dim, mystring, 0);
		
		Mat_VarSetStructFieldByName(matstruct, fieldnames[0], 0, variable); //insert in Data.name
		Mat_VarWrite(mfile.mat, matstruct, MAT_COMPRESSION_NONE);
		Mat_VarFree(matstruct);
	}
}

void DemoMatrix() {
	String filename = "matrixData.mat";
	
	MatFile mfile0;
	
	VERIFY(mfile0.OpenCreate(filename, MAT_FT_MAT73));
	
	MultiDimMatrix<double> A0(2, 3, 1);
	for (int i = 0; i < A0.GetAxisDim(0); ++i)
		for (int j = 0; j < A0.GetAxisDim(1); ++j)
			for (int k = 0; k < A0.GetAxisDim(2); ++k)
				A0(i, j, k) = i + 2*j + 3*k;
	
	mfile0.Set<double>("A", A0);

	MultiDimMatrix<std::complex<double>> Fex0(1, 3, 2);
	for (int i = 0; i < Fex0.GetAxisDim(0); ++i)
		for (int j = 0; j < Fex0.GetAxisDim(1); ++j)
			for (int k = 0; k < Fex0.GetAxisDim(2); ++k)
				Fex0(i, j, k) = std::complex<double>(i + 2*k, j + 2*k);
	
	mfile0.Set<double>("Fex", Fex0);
	
	mfile0.Set("Legend", "Hello world!");
	mfile0.SetCell("LegendList", {"First", "Second", "Third"});
	
	mfile0.Close();
	
	MatFile mfile1;
	
	VERIFY(mfile1.OpenRead(filename));	

	MultiDimMatrix<double> A1;
	mfile1.Get("A", A1);

	MultiDimMatrix<std::complex<double>> Fex1;
	mfile1.Get("Fex", Fex1);
	
	double mxA = 0;
	for (int i = 0; i < A0.size(); ++i)
		mxA = max(mxA, abs(A0.begin()[i] - A1.begin()[i]));
	VERIFY(mxA < 0.00001);

	double mxFex = 0;
	for (int i = 0; i < Fex0.size(); ++i)
		mxFex = max(mxFex, abs(Fex0.begin()[i] - Fex1.begin()[i]));
	VERIFY(mxFex < 0.00001);
	
	String legend = mfile1.Get<String>("Legend");
	VERIFY(legend == "Hello world!");
	
	Vector<String> cell;
	mfile1.GetCell("LegendList", cell);
	VERIFY(cell[0] == "First");
	VERIFY(cell[1] == "Second");
	VERIFY(cell[2] == "Third");
}

CONSOLE_APP_MAIN
{
	try {
		DemoStructure();
		DemoMatrix();
		Cout() << "\nProgram ended OK";
	} catch(Exc e) {
		Cout() << "\nError: " << e;
		SetExitCode(1);
	}

	#ifdef flagDEBUG
	Cout() << "\nPress enter to end";
	ReadStdIn();
	#endif
}
