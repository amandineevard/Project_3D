#pragma once
#include "projdyn_types.h"

namespace ProjDyn {

    struct Edge {
        Index v1, v2;
        Index t1, t2, vOtherT1, vOtherT2;
    };

    typedef std::vector<Edge> VertexStar;

    static std::vector<VertexStar> makeVertexStars(int numVertices, const Triangles& triangles) {
        std::vector<VertexStar> vertexStars;
        vertexStars.resize(numVertices);

        for (Index t = 0; t < triangles.rows(); t++) {
            for (Index v = 0; v < 3; v++) {
                Index vInd = triangles(t, v);
                for (Index ov = 0; ov < 3; ov++) {
                    if (v == ov) continue;
                    Index nbVInd = triangles(t, ov);
                    bool found = false;
                    for (Edge& checkEdge : vertexStars.at(vInd)) {
                        if (nbVInd == checkEdge.v2) {
                            checkEdge.t2 = t;
                            checkEdge.vOtherT2 = triangles(t, (Index)3 - (v + ov));
                            found = true;
                        }
                    }
                    if (!found) {
                        Edge e;
                        e.v1 = vInd;
                        e.v2 = nbVInd;
                        e.t1 = t;
                        e.vOtherT1 = triangles(t, (Index)3 - (v + ov));
                        e.t2 = -1;
                        e.vOtherT2 = -1;
                        vertexStars.at(vInd).push_back(e);
                    }
                }
            }
        }

        return vertexStars;
    }

    static std::vector<Index> getStarIndices(const VertexStar& star) {
        std::vector<ProjDyn::Index> starIndices;
        if (star.size() <= 0) return starIndices;
        starIndices.push_back(star[0].v1);
        for (int j = 0; j < star.size(); j++) starIndices.push_back(star[j].v2);
        return starIndices;
    }

	static Vector3 getTriangleNormal(const Triangles& triangles, const Positions& positions, bool dontNormalize = false)
	{
		Vector3 normal;
		normal.setZero();
		for (int i = 0; i < triangles.rows(); i++) {
			normal += (positions.row(triangles(i, 1)) - positions.row(triangles(i, 0))).cross(positions.row(triangles(i, 2)) - positions.row(triangles(i, 0))).normalized();// *triangleArea(triangles.row(i), positions);
		}
		if (!dontNormalize) normal.normalize();
		return normal;
	}


	static Scalar triangleArea(const Positions& positions, Triangle tri)
	{
		Scalar l1 = (positions.row(tri(0, 0)) - positions.row(tri(0, 1))).norm();
		Scalar l2 = (positions.row(tri(0, 1)) - positions.row(tri(0, 2))).norm();
		Scalar l3 = (positions.row(tri(0, 2)) - positions.row(tri(0, 0))).norm();
		Scalar p = (Scalar)(1. / 2.) * (l1 + l2 + l3);
		return std::sqrt(p * (p - l1) * (p - l2) * (p - l3));
	}

	static Scalar tetrahedronVolume(const Positions& positions, Tetrahedron tet) {
		// 3d edges of tet
		Eigen::Matrix<Scalar, 3, 3> edges;
		edges.col(0) = (positions.row(tet(0, 1)) - positions.row(tet(0, 0)));
		edges.col(1) = (positions.row(tet(0, 2)) - positions.row(tet(0, 0)));
		edges.col(2) = (positions.row(tet(0, 3)) - positions.row(tet(0, 0)));

		// Weight gets multiplied by tet volume
		return edges.determinant() / 6.0f;
	}

	static Scalar clamp(const Scalar& x, const Scalar& lower, const Scalar& upper) {
		return std::fmin(std::fmax(x, lower), upper);
	}

	constexpr Scalar PROJDYN_MIN_MASS = 1e-9;

	static Vector vertexMasses(const Positions& positions, const Triangles& tris) {
		Vector vMasses(positions.rows());
		vMasses.fill(0);
		int numTris = tris.rows();
		for (int tInd = 0; tInd < numTris; tInd++) {
			Scalar curArea = triangleArea(positions, tris.row(tInd)) * (1. / 3.);

			vMasses(tris(tInd, 0), 0) += curArea;
			vMasses(tris(tInd, 1), 0) += curArea;
			vMasses(tris(tInd, 2), 0) += curArea;
		}

		for (int vInd = 0; vInd < vMasses.rows(); vInd++) {
			if (vMasses(vInd, 0) < PROJDYN_MIN_MASS) {
				vMasses(vInd, 0) = PROJDYN_MIN_MASS;
			}
		}

		return vMasses;
	}

    static std::string floatToString(float v) {
        char text[50];
        int n = sprintf(text, "%f", v);
        std::string str(text, n);
        return str;
    }

}