// DGP 2019 Project
// ShapeUp and Projective Dynamics
// Author: Christopher Brandt

// This file contains the Simulator class, which implements the
// local-global algorithm that is at the heart of ShapeUp and
// Projective Dynamics.
// This class is kept as minimal as possible, and it should be
// possible to clearly see how local and global steps are
// implemented, see the methods initializeSystem() and step(),
// which are the central methods of this class.

#pragma once

#include "projdyn_types.h"
#include "projdyn_common.h"
#include "projdyn_constraints.h"

#include <memory>
#include <iostream>

constexpr double PROJDYN_INITIAL_STIFFNESS = 100.;

namespace ProjDyn {
	typedef Eigen::SimplicialLDLT<SparseMatrix> SparseSolver;

	// The simulator class, which handles objects represented as Eigen matrices
	// and the constraints above to run a Projective Dynamics simulation.
	class Simulator {
	public:

		/*
			Creates an empty Simulator, that can be re-used for multiple simulations.
			To set it up you will need to set a mesh (setMesh()), add constraints
			(addConstraints()) and then initialize the system (initializeSystem()).
			After that, the simulation can be run by using the step() function for
			each frame.
			Updated vertex positions can then be extracted using getPositions().
		*/
		Simulator(Scalar time_step = 1. / 60.) :
			m_positions(0, 3),
			m_initial_positions(0, 3),
			m_triangles(0, 3),
			m_tetrahedrons(0, 4),
			m_mass(1.),
			m_stiffness_factor(PROJDYN_INITIAL_STIFFNESS),
			m_mass_matrix(0, 0),
			m_vertex_masses(0, 1),
			m_num_verts(0),
			m_num_tris(0),
			m_num_tets(0),
			m_time_step(time_step),
			m_lhs(0, 0),
			m_constraint_mat_t(0, 0),
			m_laplacian(0, 0),
			m_system_init(false),
			m_lhs_updated(false)
		{}

        // Based on the current vertices, triangles, tetrahedrons and constraints,
        // as well as the current time-step and stiffness-factor, the global system
        // will be set up (as well as some structures for the local steps), such
        // that the simulation can be run by using the step() function.
        // NOTE: the method setMesh() needs to be called before runnin this method,
        // otherwise, no mesh is available.
        bool initializeSystem() {
            // Check if vertices and triangles or tetrahedrons are available
            // Check if constraints are available
            if (m_num_verts == 0 || (m_num_tris == 0 && m_num_tets == 0)) {
                std::cout << "WARNING: Could not initialize: either no vertices, or no faces and no tets." << std::endl;
                return false;
            }

            // Get list of all current constraints from constraint group
            m_constraints.clear();
            for (auto cg : m_constraintGroups) {
                for (auto c : cg->constraints) {
                    // Update weight multiplier of constraint and add to list of all constraints
                    c->setWeightMultiplier(cg->weight);
                    m_constraints.push_back(c);
                }
            }

            if (m_constraints.size() == 0) {
                std::cout << "WARNING: Could not initialize: no constraints." << std::endl;
                return false;
            }

            // Collect entries of constraint matrix
            std::vector<Triplet> con_triplets, con_triplets_t;
            Index row_ind = 0, row_ind2 = 0;;
            for (auto c : m_constraints) {
                // The following function adds the row sqrt(w_i) A_i S_i
                // (using the notation of the paper) to the matrix
                c->addConstraint(con_triplets, row_ind2, true, false);

                // The following function adds the row w_i S_i^T A_i^T
                // (using the notation of the paper) to the matrix
                c->addConstraint(con_triplets_t, row_ind, false, true);
            }

            // Construct the "laplacian" matrix from the triplets,
            // i.e. the matrix Sum_i w_i S_i^T A_i^T A_i S_i, which
            // appears on the left hand side of the system of the global step
            SparseMatrixRM temp1(row_ind, m_num_verts);
            temp1.setFromTriplets(con_triplets.begin(), con_triplets.end());
            SparseMatrix temp2 = temp1.transpose();
            m_laplacian = temp2 * temp1;

            // Also construct the matrix Sum_i w_i S_i^T A_i^T 
            // which appears on the right hand side of the global step
            m_constraint_mat_t.resize(m_num_verts, row_ind);
            m_constraint_mat_t.setFromTriplets(con_triplets_t.begin(), con_triplets_t.end());

            // Initialize internal quantities for simulation
            m_momentum.resize(m_num_verts, 3);
            m_momentum.setZero();
            m_velocities.resize(m_num_verts, 3);
            m_velocities.setZero();
            m_ext_forces.resize(m_num_verts, 3);
            m_ext_forces.setZero();
            m_ext_forces.col(1).setConstant(-m_gravity);
            m_constraint_projections.resize(m_constraint_mat_t.cols(), 3);
            m_constraint_projections.setZero();
            m_rhs.resize(m_num_verts, 3);
            m_rhs.setZero();
            m_old_positions.resize(m_num_verts, 3);
            m_old_positions.setZero();

            // With this, the system is initialized
            m_system_init = true;

            // Compute and factorize the lhs matrix
            return recomputeLHS();
        }

        // Performs a step by applying the specified amount of local global iterations.
        // The system needs to be initialized and the lhs matrix up-to-date.
        // After these steps, the updated positions can be extracted via getPositions()
        // NOTE: for ShapeUp style shape editing, m_dynamicMode is set to false, and all
        // sections that require it to be true can be ignored.
        bool step(int num_iterations) {
            // The system needs to be initialized before running the local-global alg.
            if (!m_system_init || !m_lhs_updated) return false;

            if (m_dynamicMode) {
                // Store old positions to later compute the velocities
                m_old_positions = m_positions;
                // Compute the momentum term, which will be fixed for all iterations
                m_momentum = m_positions + m_time_step * m_velocities + m_time_step * m_time_step * m_ext_forces;
                m_positions = m_momentum;
            }
            else {
                m_momentum = m_positions;
            }

            // If vertices are being grabbed by the mouse, we enforce this here 
            // (this is only used in the simulation, for ShapeUp, proper position
            // constraint groups are used!)
            if (m_hasGrab && m_grabVerts.size() == m_grabPos.size()) {
                for (int j = 0; j < m_grabVerts.size(); j++) {
                    Index ind = m_grabVerts[j];
                    Eigen::Vector3f pos = m_grabPos[j];
                    for (int d = 0; d < 3; d++) m_momentum(ind, d) = pos(d);
                }
            }

            // The local-global algorithm:
            for (int it = 0; it < num_iterations; it++) {

                // --------------------------------------------------
                // Local step: compute constraint projections and put them in a vector
                // --------------------------------------------------
#pragma omp parallel for
                for (int j = 0; j < m_constraints.size(); j++) {
                    m_constraints[j]->project(m_positions, m_constraint_projections);
                }

                // --------------------------------------------------
                // Global step: construct rhs and solve linear system
                // --------------------------------------------------
                // Compute the r.h.s of the linear system
                m_rhs = m_stiffness_factor * m_constraint_mat_t * m_constraint_projections;
                if (m_dynamicMode) {
                    m_rhs += ((Scalar)1 / (m_time_step * m_time_step)) * m_mass_matrix * m_momentum;
                }
                // Solve the three linear systems cooresponding to the x, y, and z coordinate
#pragma omp parallel for
                for (int coord = 0; coord < 3; coord++) {
                    m_positions.col(coord) = m_solver.solve(m_rhs.col(coord));
                }
            }

            if (m_dynamicMode) {
                // Lastly, update the velocities as finite differences between the old and new positions
                m_velocities = (1. / m_time_step) * (m_positions - m_old_positions);
            }

            return true;
        }

		// Provide a n by 3 scalar matrix containing x, y, z positions of each vertex per row
		// and a m by 3 index matrix containing the indices of the vertices in each triangle
		void setMesh(Positions& pos, Triangles& tris) {
            clearConstraints();

			m_positions = pos;
			m_initial_positions = pos;
			m_triangles = tris;

			// Set the floor below the mesh at a distance corresponding to twice the height
			// of the bounding box
			m_floorHeight = m_positions.col(1).minCoeff() - (m_positions.col(1).maxCoeff() - m_positions.col(1).minCoeff()) * 2;
			meshChanged();
		}

		// Provide a n by 3 scalar matrix containing x, y, z positions of each vertex per row
		// and a m by 3 index matrix containing the indices of the vertices in each OUTER (boundary) triangle
		// and a l by 4 index matrix containing the indices of the vertices all tetrahedrons
		void setMesh(Positions& pos, Triangles& tris, Tetrahedrons& tets) {
			m_tetrahedrons = tets;
			setMesh(pos, tris);
		}

		// Reset vertex positions to their initial positions, as given when setMesh() was called.
		void resetPositions() {
			m_positions = m_initial_positions;
			m_velocities.setZero(m_positions.rows(), 3);
		}

		const Triangles& getTriangles() const {
			return m_triangles;
		}

		const Tetrahedrons& getTetrahedrons() const {
			return m_tetrahedrons;
		}

		// Positions of the deformed mesh
		const Positions& getPositions() const {
			return m_positions;
		}

		// Positions of the rest-state mesh
		const Positions& getInitialPositions() const {
			return m_initial_positions;
		}

		// Returns the number of outer vertices.
		// Note that the vertices with index 0, ..., numOuterVerts-1 are the vertices on
		// the surface triangles of the mesh, while the rest are inner vertices.
		Index getNumOuterVerts() const {
			return m_num_outer_verts;
		}

		// Returns the total number of vertices (inner and outer vertices)
		Index getNumVerts() const {
			return m_num_verts;
		}

		// Adds constraints to the simulation.
		// The system will have to be re-initialized to contain the changed constraints.
		void addConstraints(const std::vector<ConstraintPtr>& newCons) {
            ConstraintGroupPtr newGroup = std::make_shared<ConstraintGroup>("Group " + m_globalIDCounter, newCons, 1.);
            m_globalIDCounter++;
            m_constraintGroups.push_back(newGroup);
			m_system_init = false;
		}

        void addConstraints(ConstraintGroupPtr newCons) {
            m_constraintGroups.push_back(newCons);
            m_system_init = false;
        }

        void removeConstraints(ConstraintGroupPtr constraints) {
            m_constraintGroups.erase(std::remove(m_constraintGroups.begin(), m_constraintGroups.end(), constraints), m_constraintGroups.end());
        }

        void removeConstraints(std::string name) {
            auto remove_it = std::remove_if(m_constraintGroups.begin(), m_constraintGroups.end(), [name](const ConstraintGroupPtr& v) { return (v->name == name); });
            m_constraintGroups.erase(remove_it, m_constraintGroups.end());
        }

        // Adds single constraint to the simulation.
        // The system will have to be re-initialized to contain the changed constraints.
        void addConstraint(const ConstraintPtr& con) {
            ConstraintGroupPtr newGroup = std::make_shared<ConstraintGroup>("Group " + m_globalIDCounter, std::vector<ConstraintPtr>({ con }), 1.);
            m_globalIDCounter++;
            m_constraintGroups.push_back(newGroup);
            m_system_init = false;
        }

		// Removes all constraints from the simulation.
		void clearConstraints() {
			m_constraintGroups.clear();
			m_system_init = false;
		}

		const std::vector<ConstraintGroupPtr>& getConstraintGroups() {
			return m_constraintGroups;
		}

		// Only if this returns true will the step() method perform a simulation
		// taking into account the current mesh and constraints.
		// Otherwise, initializeSystem() needs to be called first.
		bool isInitialized() const {
			return (m_system_init && m_lhs_updated);
		}

		// Set gravitational acceleration.
		// Note that external forces are reset to 0 when the system is initialized.
		void setGravity(Scalar g) {
            m_gravity = g;
		}

        // Switches between dynamic (i.e. simulation) mode and static mode (i.e. 
        // constraint based shape optimization)
        void setDynamic(bool dynamic) {
            m_dynamicMode = dynamic;
        }

		// Builds and factorizes the linear system.
		// Needs to be re-run everytime the stiffness-factor or time-step changes,
		// and is automatically run each time initializeSystem is called.
		bool recomputeLHS() {
			if (!m_system_init) return false;
			m_lhs = m_stiffness_factor * m_laplacian;
            if (m_dynamicMode) {
                m_lhs += ((Scalar)1. / (m_time_step * m_time_step)) * m_mass_matrix;
            }
			m_solver.compute(m_lhs);
			m_lhs_updated = (m_solver.info() == Eigen::Success);
			if (!m_lhs_updated) {
				std::cout << "WARNING: could not initialize system, LHS matrix not pos.def.? (Possibly not all vertices take part in constraints?)" << std::endl;
            }
            else {
                std::cout << "Successfully recomputed and refactorized LHS of global step" << std::endl;
            }

			return m_lhs_updated;
		}

		// Change the time-step.
		// Since this changes the system matrix, it sets the lhs as not updated
		void setTimeStep(Scalar time_step) {
			m_time_step = time_step;
			m_lhs_updated = false;
		}

		// Change the stiffness-factor
		// Since this changes the system matrix, it sets the lhs as not updated
		void setStiffnessFactor(Scalar stiffness_factor) {
			m_stiffness_factor = stiffness_factor;
			m_lhs_updated = false;
		}

		void setGrab(const std::vector<Index>& grabVerts, const std::vector<Eigen::Vector3f> grabPos) {
			m_grabVerts = grabVerts;
			m_grabPos = grabPos;
			m_hasGrab = true;
		}

		void releaseGrab() {
			m_hasGrab = false;
		}

		void addFloorConstraints(Scalar weightMultiplier, Scalar forceFactor = 1.) {
            if (m_num_verts <= 0) return;
			Vector voronoiAreas = vertexMasses(getInitialPositions(), getTriangles());
            std::vector<ConstraintPtr> floorCons;
			for (Index v = 0; v < m_num_verts; v++) {
				floorCons.push_back(std::make_shared<FloorConstraint>(v, voronoiAreas(v) * weightMultiplier, m_floorHeight, forceFactor));
			}
            addConstraints(std::make_shared<ConstraintGroup>("Floor", floorCons, 1));
			m_system_init = false;
		}

        Scalar getFloorHeight() const {
            return m_floorHeight;
        }

	protected:
		// Mesh faces, vertices and tetrahedrons
		Positions m_positions, m_initial_positions;
		Triangles m_triangles;
		Tetrahedrons m_tetrahedrons;
		Index m_num_verts = 0;
		Index m_num_outer_verts = 0;
		Index m_num_tris = 0;
		Index m_num_tets = 0;

		// External forces per vertex
		// Note that those are pre-multiplied by the inverse mass, such that
		// the momentum term can be computed cheaper.
		Positions m_ext_forces;

		// Internal quantities during simulation
		Positions m_velocities, m_momentum, m_old_positions;
		Positions m_constraint_projections;
		Positions m_rhs;

		// Mass per unit area (or unit volume, if tets are provided)
		Scalar m_mass = 1.;

		// Vertex masses and mass matrix
		SparseMatrix m_mass_matrix;
		Vector m_vertex_masses;

		// Stiffness factor, applied to *all* constraints
		Scalar m_stiffness_factor = 1.;

		// List of constraint groups
		std::vector<ConstraintGroupPtr> m_constraintGroups;

        // List of actual constraints, filled from constrain groups list when calling initializeSystem
        std::vector<ConstraintPtr> m_constraints;

		// Transposed constraint matrix sum w_i S_i A_i.
		// Defined as row-major, since constraints are listed row-by-row
		SparseMatrix m_constraint_mat_t;

		// The "laplacian" sum w_i A_i^T S_i^T S_i A_i
		SparseMatrix m_laplacian;

		// The lhs matrix of the global step, i.e. (1/h^2) M + s * A^T A,
		// where h is the time-step and s is the stiffness factor
		SparseMatrix m_lhs;

		// The linear solver for the global step, contains the factorization, if m_lhs_updated == true
		SparseSolver m_solver;

		// Will be set to false when time-step or stiffness factor change, which requires
		// a recomputation of the lhs
		bool m_lhs_updated = false;

		// If false, the system needs to be set up and factorized before being able to simulate
		bool m_system_init = false;

		// Time step of the simulation
		Scalar m_time_step = 1. / 60.;

		// If m_hasGrab is true, vertices with the indices in m_grabVerts will be forced
		// to the positions set in m_grabPos
		bool m_hasGrab = false;
		std::vector<Index> m_grabVerts;
		std::vector<Eigen::Vector3f> m_grabPos;

		// y-Coordinate of the floor plane, used for preliminary collisions
		Scalar m_floorHeight = 0;

        // Toggles "Dynamic mode", i.e. if set to true, this object is used to simulate
        // an object under external forces and momentum, if set to false, momentum and
        // external forces are ignored and this object is used to solve shape optimization
        // problems
        bool m_dynamicMode = true;

        // Gravity is added as an external force acting on the y-component of the system
        Scalar m_gravity = 0;

        // Simply to give a consecutive numbering of constraint groups
        Index m_globalIDCounter = 0;

		// Reset constraints, update mass matrix, set system status as not initialized
		void meshChanged() {
			m_num_verts = m_positions.rows();
			m_num_tris = m_triangles.rows();
			m_num_tets = m_tetrahedrons.rows();
			m_num_outer_verts = m_triangles.maxCoeff() + 1;
			m_constraints.clear();
			computeMassMatrix();
			m_system_init = false;
		}

		void computeMassMatrix() {
			m_vertex_masses.setZero(m_num_verts);
			// Compute masses
			if (m_num_tets == 0) {
				// Surface mesh: compute vertex masses from face areas
				for (int i = 0; i < m_num_tris; i++) {
					Scalar cur_mass = triangleArea(m_positions, m_triangles.row(i));
					for (int j = 0; j < 3; j++) {
						m_vertex_masses(m_triangles(i, j)) += cur_mass / 3;
					}
				}
			}
			else {
				// Volume mesh: compute vertex masses from face areas
				for (int i = 0; i < m_num_tets; i++) {
					Scalar cur_mass = tetrahedronVolume(m_positions, m_tetrahedrons.row(i));
					for (int j = 0; j < 4; j++) {
						m_vertex_masses(m_tetrahedrons(i, j)) += cur_mass / 4;
					}
				}
			}
			m_mass_matrix.resize(m_num_verts, m_num_verts);
			m_mass_matrix.setIdentity();
			for (int i = 0; i < m_num_verts; i++) {
				m_mass_matrix.coeffRef(i, i) = m_vertex_masses(i);
			}
		}
	};

}
