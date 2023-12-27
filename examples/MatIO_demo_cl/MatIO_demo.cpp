// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#include <Core/Core.h>
#include "data.brc"

using namespace Upp;

#include <MatIO/matio.h>

CONSOLE_APP_MAIN
{
	MatFile mat;
	
	String fileName = "temp_file_2.mat";
	
	if (!FileExists(fileName)) 
		VERIFY(SaveFile(fileName, String(dataMat, dataMat_length))); // Error: Impossible to save .mat file
	
	VERIFY(mat.OpenWrite(fileName));	 
		
	Cout() << "\nVersion " << mat.GetVersionName();
	
	Cout() << "\nVariables list:";
	for (int i = 0; i < mat.GetVarCount(); ++i) {
		MatVar var = mat.GetVar(mat.GetVarName(i));
		Cout() << Format("\n- %s <%s>(%d, %d)", var.GetName(), var.GetTypeString(), var.GetDimCount(0), var.GetDimCount(1));
	}
	
	MatVar est = mat.GetVar("Structure");
	if (est.IsLoaded()) {
		Cout() << "\nStructure " << est.GetTypeString();
		int numf = est.GetFieldCount();
		for (int i = 0; i < numf; ++i) {
			String name = est.GetFieldName(i);
			MatVar var = est.GetVar(name);
			Cout() << Format("\n- %s <%s>(%d, %d)", name, var.GetTypeString(), var.GetDimCount(0), var.GetDimCount(1));
		}
		Cout() << "\nStructure.Est2.One: " << mat.VarRead<int>(mat.GetVar("Structure").GetVar("Est2").GetVar("One"));
		Cout() << "\nStructure.Est2.Two: " << mat.VarRead<String>(mat.GetVar("Structure").GetVar("Est2").GetVar("Two"));
		
		mat.VarDelete("Estructura");
	}
	
	Vector<String> elements;
	elements << "First" << "Est2" << "Tercero";
	MatVar estructura("Estructura", 1, 1, elements);
	estructura.VarWriteStruct("First", "The first");

	Vector<String> elements2;
	elements2 << "One" << "Two";
	MatVar est2("Est2", 1, 1, elements2);
	est2.VarWriteStruct("Two", "The second 22");
	est2.VarWriteStruct("One", 211111);
	
	estructura.VarWriteStruct("Est2", est2);
	
	mat.VarWrite(estructura);
	

/*	MatVar opt = mat.GetVar("Options");
	Cout() << "\nOptions " << opt.GetTypeString();
	int numf = opt.GetFieldCount();
	for (int i = 0; i < numf; ++i) {
		String name = opt.GetFieldName(i);
		MatVar var = opt.GetVar(name);
		Cout() << Format("\n- %s <%s>(%d, %d)", name, var.GetTypeString(), var.GetDimCount(0), var.GetDimCount(1));
	}
	
	Cout() << "\nFreqRangeChoice: " << mat.VarReadString(mat.GetVar("Options").GetVar("FreqRangeChoice"));
	Cout() << "\nFreqChoice: "		<< mat.VarReadString(mat.GetVar("Options").GetVar("FreqChoice"));
	Cout() << "\nMode: " 			<< mat.VarRead<double> (mat.GetVar("Options").GetVar("Mode"));
	Cout() << "\nMod: " 			<< mat.VarRead<double> (mat.GetVar("Options").GetVar("Mod"));


	MatVar optim = opt.GetVar("Optim");
	int numOpt = optim.GetFieldCount();
	for (int i = 0; i < numOpt; ++i) {
		String name = optim.GetFieldName(i);
		MatVar var = optim.GetVar(name);
		Cout() << Format("\n- %s <%s>(%d, %d)", name, var.GetTypeString(), var.GetDimCount(0), var.GetDimCount(1));
	}
	
	Cout() << "\nInit.Cond: " 	  << mat.VarRead<double>(mat.GetVar("Options").GetVar("Optim").GetVar("InitCond"));
	Cout() << "\nInit.maxEval: "  << mat.VarRead<double>(mat.GetVar("Options").GetVar("Optim").GetVar("maxEval"));
	
	Cout() << "\nInit.Tol: " 	  << mat.VarRead<double>(mat.GetVar("Options").GetVar("Optim").GetVar("Tol"));
	Cout() << "\nInit.StepTol: "  << mat.VarRead<double>(mat.GetVar("Options").GetVar("Optim").GetVar("StepTol"));
	Cout() << "\nInit.ThresRel: " << mat.VarRead<double>(mat.GetVar("Options").GetVar("Optim").GetVar("ThresRel"));
	Cout() << "\nInit.ThresAbs: " << mat.VarRead<double>(mat.GetVar("Options").GetVar("Optim").GetVar("ThresAbs"));

	mat.VarDelete("Options");
*/	
	/*
	//Create a 1 x 1 struct 'Data' with fields: name, unit, value
	char *structname = "Data";
	const char *fieldnames[3] = { "name","unit","value" };
	size_t structdim[2] = { 1, 1 }; // create 1 x 1 struct
	matvar_t* matstruct = Mat_VarCreateStruct(structname, 2, structdim, fieldnames, 3); 

	char* mystring = "Speed";
	size_t dim[2] = { 1, 5 }; //string dimension
	matvar_t *variable = Mat_VarCreate(fieldnames[0], MAT_C_CHAR, MAT_T_UTF8, 2, dim, mystring, 0);
	
	Mat_VarSetStructFieldByName(matstruct, fieldnames[0], 0, variable); //insert in Data.name
	Mat_VarWrite(mat.mat, matstruct, MAT_COMPRESSION_NONE);
	Mat_VarFree(matstruct);
*/

	
	MatMatrix<double> aa = mat.VarReadMat<double>("A");		
	Cout() << "\nA";
	aa.Print();
	
	if (mat.VarExists("A_b")) {
		MatMatrix<double> ab = mat.VarReadMat<double>("A_b");
		Cout() << "\nA_b";
		ab.Print();
	}
	
	aa(3, 0) = 5.8;
	
	Cout() << "\nA_b";
	aa.Print();

	VERIFY(mat.VarWrite("A_b", aa));
	
	MatMatrix<double> w = mat.VarReadMat<double>("w");		
	Cout() << "\nw";
	w.Print();
	
	Cout() << "\nProgram ended";
	#ifdef flagDEBUG
	ReadStdIn();
	#endif
}
