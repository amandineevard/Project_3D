#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include "projdyn_types.h"
#include "tetgen.h"

namespace ProjDyn {

	// Small interface to tetgen that allows for meshes represented by Eigen matrices to be tetrahedralized, keeping their
	// outer vertices and faces the same, while adding inner points and tetrahedra.
	// The matrices verts and tris are input, while newVertsOut and tetsOut are output.
	// q is a quality parameter that defines the maximal ratios of edge lengths within a tet.
	static bool tetrahedralizeMesh(Positions& verts, Triangles& tris, Positions& newVertsOut, Tetrahedrons& tetsOut, float q = 1.1) {
		std::ofstream polyFile;
		std::string fileName = "temp.poly";
		if (tris.rows() == 0) {
			fileName = "temp.node";
		}
		polyFile.open(fileName);

		polyFile << verts.rows() << " 3 0 0\n";
		for (Index i = 0; i < verts.rows(); i++) {
			polyFile << (i + 1) << " " << verts(i, 0) << " " << verts(i, 1) << " " << verts(i, 2) << "\n";
		}

		if (tris.rows() > 0) {
			polyFile << tris.rows() << " 0\n";
			for (Index i = 0; i < tris.rows(); i++) {
				polyFile << "1\n3 " << (tris(i, 2) + 1) << " " << (tris(i, 1) + 1) << " " << (tris(i, 0) + 1) << "\n";
			}
			polyFile << "0\n0\n";
		}

		polyFile.close();

		std::ostringstream argsStr;
		argsStr << "-pq" << q;
		std::string args = argsStr.str();

		if (tris.rows() == 0) args = "";

		std::cout << "Arguments used for tetgen: " << args << std::endl;

		std::string nothing = " ";

		// Turn arguments and filename into input for the tetgen tetrahedralize() method
		char* argsC = new char[args.length() + 1];
		strcpy(argsC, args.c_str());
		char* fileC = new char[fileName.length() + 1];
		strcpy(fileC, fileName.c_str());
		char* nothingC = new char[nothing.length() + 1];
		strcpy(nothingC, nothing.c_str());
		char** argv = new char* [3]{ nothingC, argsC, fileC };
		tetgenbehavior b;
		tetgenio in, addin, bgmin;
		if (!b.parse_commandline(3, argv)) {
			terminatetetgen(NULL, 10);
			return false;
		}
		// Read input files.
		if (b.refine) { // -r
			if (!in.load_tetmesh(b.infilename, (int)b.object)) {
				terminatetetgen(NULL, 10);
				return false;
			}
		}
		else { // -p
			if (!in.load_plc(b.infilename, (int)b.object)) {
				terminatetetgen(NULL, 10);
				return false;
			}
		}
		if (b.insertaddpoints) { // -i
								 // Try to read a .a.node file.
			addin.load_node(b.addinfilename);
		}
		if (b.metric) { // -m
						// Try to read a background mesh in files .b.node, .b.ele.
			bgmin.load_tetmesh(b.bgmeshfilename, (int)b.object);
		}

		// Start tetrahedralization using tetgen
		tetrahedralize(&b, &in, NULL, &addin, &bgmin);

		// Create new vertices matrix
		newVertsOut = verts;
		std::ifstream nodeFile;
		nodeFile.open("temp.1.node");
		if (nodeFile.is_open()) {
			std::cout << "Getting new vertices from tet-mesh..." << std::endl;

			std::string line;
			std::getline(nodeFile, line);
			std::istringstream iss(line);
			Index numVertsNew;
			iss >> numVertsNew;
			if (numVertsNew < verts.rows()) {
				std::cout << "Error: New node file contains less vertices than the mesh!" << std::endl;
				return false;
			}
			else if (numVertsNew > verts.rows()) {
				newVertsOut.conservativeResize(numVertsNew, 3);
				while (std::getline(nodeFile, line))
				{
					std::istringstream iss(line);
					Index ind;
					float x, y, z;
					if (!(iss >> ind >> x >> y >> z)) {
						std::cout << "Error: Node file has wrong format, couldn't read vertex " << ind << "!" << std::endl;
						return false;
					}

					ind--; // Vertices in tetgen files start from index 1

					if (ind >= verts.rows()) {
						newVertsOut.row(ind) = Vector3(x, y, z);
					}

					if (ind >= numVertsNew - 1) {
						break;
					}
				}
			}
			nodeFile.close();
		}
		else {
			std::cout << "Error: Couldn't find .node (new vertices) file " << std::endl;
			return false;
		}

		// Create new tets matrix
		std::ifstream tetsFile;
		tetsFile.open("temp.1.ele");
		if (tetsFile.is_open()) {
			std::cout << "Getting new tets from tet-mesh..." << std::endl;

			std::string line;
			std::getline(tetsFile, line);
			std::istringstream iss(line);
			Index numTets;
			iss >> numTets;

			tetsOut.resize(numTets, 4);
			Index tetInd = 0, lineInd = 0, skipped = 0;
			while (std::getline(tetsFile, line))
			{
				std::istringstream iss(line);
				Index tetIndTG, ind1, ind2, ind3, ind4;
				if (!(iss >> tetIndTG >> ind1 >> ind2 >> ind3 >> ind4)) {
					std::cout << "Error: Tet file has wrong format, couldn't read tet " << tetInd << "!" << std::endl;
					return false;
				}

				if (ind1 != ind2 && ind1 != ind3 && ind1 != ind4 && ind2 != ind3 && ind2 != ind4 && ind3 != ind4) {
					// Vertices in tetgen files start from index 1!!
					tetsOut(tetInd, 0) = ind1 - 1;
					tetsOut(tetInd, 1) = ind2 - 1;
					tetsOut(tetInd, 2) = ind3 - 1;
					tetsOut(tetInd, 3) = ind4 - 1;
					tetInd++;
				}
				else {
					skipped++;
				}

				if (lineInd >= numTets - 1) {
					break;
				}

				lineInd++;
			}

			tetsOut.conservativeResize(tetInd, 4);

			if (skipped > 0) {
				std::cout << "Skipped " << skipped << " unwanted tetrahedra." << std::endl;
			}

			tetsFile.close();
		}
		else {
			std::cout << "Error: Couldn't find .ele (tet) file " << std::endl;
			return false;
		}

		return true;
	}

}