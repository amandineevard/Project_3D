// DGP 2019 Project
// ShapeUp and Projective Dynamics
// Author: Christopher Brandt

// An API that provides simple methods to handle the interaction between
// the viewer and the ProjectiveDynamics simulation

#pragma once

#include "projdyn_types.h"
#include <thread>
#include "projdyn.h"
#include "projdyn_tetgen.h"
#include <nanogui/slider.h>
#include <nanogui/checkbox.h>
#include "viewer.h"
#include "projdyn_widgets.h"

class ProjDynAPI {
public:
    // Some default values
    const bool UPDATE_NORMALS = true;
    const int  NUM_ITS_INITIAL = 10;
    const int  FPS = 60;
    const bool DYNAMIC_MODE = true;

    ProjDynAPI(Viewer* viewer) {
        m_numIterations = NUM_ITS_INITIAL;
        setDynamic(DYNAMIC_MODE);
        m_viewer = viewer;
    }

    // Changes the behavior of the simulator:
    // If it is set to non-dynamic, it behaves like a ShapeUp style shape defomration
    // framework.
    // If it is set to dynamics, it behaves like a Projective Dynamics simulation.
    void setDynamic(bool dynamic) {
        m_simulator.setDynamic(dynamic);
    }

    // Create the GUI for controlling a Projective Dynamics simulation:
    // starting, pausing and resetting a simulation, as well as adding
    // constraints, defining the number of iterations, etc.
    void initSimulationGUI(bool addConstraintGUI = false) {
        Window* pd_win = new Window(m_viewer, "Simulation Controls");
        pd_win->setPosition(Vector2i(15, 230));
        pd_win->setLayout(new GroupLayout());

        Widget* panel = new Widget(pd_win);
        panel->setLayout(new BoxLayout(nanogui::Orientation::Vertical, nanogui::Alignment::Middle, 0, 10));

        m_startButton = new Button(panel, "Run Simulation");
        m_startButton->setFlags(Button::RadioButton);
        m_startButton->setPushed(false);
        m_startButton->setCallback([this]() {
            start();
        });

        m_stopButton = new Button(panel, "Stop Simulation");
        m_stopButton->setFlags(Button::RadioButton);
        m_stopButton->setPushed(true);
        m_stopButton->setCallback([this]() {
            stop();
        });

        CheckBox* updateNormalsCB = new CheckBox(pd_win, "Update Normals");
        updateNormalsCB->setChecked(m_updateNormals);
        updateNormalsCB->setCallback([this](bool state) {
            m_updateNormals = state;
        });

        Button* reset_b = new Button(pd_win, "Reset Positions");
        reset_b->setCallback([this]() {
            bool was_active = m_simActive;
            stop();
            m_simulator.resetPositions();
            if (was_active) {
                start();
            }
            else {
                uploadPositions();
            }
        });


        Button* addtets_b = new Button(pd_win, "Tetrahedralize");
        addtets_b->setCallback([this]() {
            setMesh(true);
        });

        PopupButton* popupBtnTemp = new PopupButton(pd_win, "Temperature model");
        Popup* popupTemp = popupBtnTemp->popup();
        popupTemp->setLayout(new GroupLayout());

        Button* d = new Button(popupTemp, "Constant");
        d->setCallback([this, popupBtnTemp]() {
            bool was_active = m_simActive;
            stop();
            m_simulator.setTemperatureModel(uniform);
            m_simulator.setTempCoefUniform(0.0);
            updateConstraintsGUI();
			updateTemperatureGUI();
            if (was_active) {
                start();
            }
            popupBtnTemp->setPushed(false);
        });
        d = new Button(popupTemp, "Linear");
        d->setCallback([this, popupBtnTemp]() {
            bool was_active = m_simActive;
            stop();
            m_simulator.setTemperatureModel(linear);
            m_simulator.setTempCoefLinear(0.0);
            updateConstraintsGUI();
			updateTemperatureGUI();
            if (was_active) {
                start();
            }
            popupBtnTemp->setPushed(false);
        });
        d = new Button(popupTemp, "Diffusion");
        d->setCallback([this, popupBtnTemp]() {
            bool was_active = m_simActive;
            stop();
            m_simulator.setTemperatureModel(diffusion);
            m_simulator.setTempCoefDiffusion(0.0);
            updateConstraintsGUI();
			updateTemperatureGUI();
            if (was_active) {
                start();
            }
            popupBtnTemp->setPushed(false);
        });


        PopupButton* popupBtn = new PopupButton(pd_win, "Add constraints", ENTYPO_ICON_LINK);
        Popup* popup = popupBtn->popup();
        popup->setLayout(new GroupLayout());

        Button* b = new Button(popup, "Floor");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            m_simulator.addFloorConstraints(10., 3.);
            m_viewer->setFloorHeight(m_simulator.getFloorHeight());
            m_viewer->showFloor(true);
            updateConstraintsGUI();
			updateTemperatureGUI();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Edge Springs");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            addEdgeSpringConstraints();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Triangle Strain");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            addTriangleStrainConstraints();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Triangle Bending");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            addBendingConstraints();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Tet Strain");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            addTetStrainConstraints();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        Button* clear_b = new Button(pd_win, "Clear constraints");
        clear_b->setCallback([this]() {
            stop();
            clearConstraints();
            m_simulator.resetPositions();
            m_viewer->showFloor(false);
            uploadPositions();
            updateConstraintsGUI();
			updateTemperatureGUI();
        });

        Label* iterations_label = new Label(pd_win, "Num Loc-Glob Its: ");
        IntBox<int>* iterations_box = new IntBox<int>(pd_win, m_numIterations);
        iterations_box->setEditable(true);
        iterations_box->setCallback([this](int num_its) {
            m_numIterations = num_its;
        });

        if (addConstraintGUI) {
            initConstraintsGUI();
        }

		initTemperatureGUI();

        m_viewer->performLayout();
    }

    // Add an (empty) GUI window in which each constraint group that is added through
    // this API gets a slider which controls a weight multiplier.
    void initConstraintsGUI() {
        m_constraint_window = new Window(m_viewer, "Constraints weights");
        m_constraint_window->setPosition(Vector2i(780, 25));
        m_constraint_window->setLayout(new GroupLayout());
    }

	// Add an (empty) GUI window in which temperature
	// gets a slider which controls a weight multiplier.
	void initTemperatureGUI() {
		m_temperature_window = new Window(m_viewer, "Temperature coefficients");
		m_temperature_window->setPosition(Vector2i(780, 360));
		m_temperature_window->setLayout(new GroupLayout());
	}

    // Each time the constraints of the simulation change, this gets called
    // to create one slider (and textboxes) for each constraint group.
    void updateConstraintsGUI() {
        if (!m_constraint_window) return;

        // Clear all previous constraint controls
        while (m_constraint_window->children().size() > 0) {
            m_constraint_window->removeChild(0);
        }

        // For each constraint group add controls
        for (const auto& g : m_simulator.getConstraintGroups()) {
            new Label(m_constraint_window, g->name, "sans-bold");
            Widget* panel = new Widget(m_constraint_window);
            panel->setLayout(new BoxLayout(nanogui::Orientation::Horizontal, nanogui::Alignment::Middle, 0, 10));

            // Add a sliderand set defaults
            ConstraintSlider* slider = new ConstraintSlider(panel, m_viewer, m_simulator.getNumVerts(), g);

            // Re-initialize system and update positions once the user lets go of the slider
            slider->setFinalCallback([this, slider](float v) {
                slider->setValue(v);
                bool wasRunning = m_simActive;
                stop();
                m_simulator.initializeSystem();
                update();
                if (wasRunning) start();
            });
        }

        Button* b = new Button(m_constraint_window, "Update");
        b->setCallback([this]() {
            update();
        });

        m_viewer->performLayout();
    }

	// Each time the temperature model or the constraints of the simulation change,
	// this gets called to create one slider (and textboxes) for the
	// temperature model and for each temperature dependent constraint group.
	void updateTemperatureGUI() {
		if (!m_temperature_window) return;

		// Clear all previous temperature controls
		while (m_temperature_window->children().size() > 0) {
			m_temperature_window->removeChild(0);
		}

		// Add temperature model coefficient slider:
		TemperatureModel tm = m_simulator.getTemperatureModel();
		std::string name = "";
		if (tm == uniform){
			name = "Uniform temperature";
		}
		else if (tm == linear){
			name = "Floor temperature";
		}
		else if (tm == diffusion){
			name = "Diffusion coefficient";
		}
		new Label(m_temperature_window, name, "sans-bold");
		Widget* panel = new Widget(m_temperature_window);
		panel->setLayout(new BoxLayout(nanogui::Orientation::Horizontal, nanogui::Alignment::Middle, 0, 10));
		TemperatureSlider* slider = new TemperatureSlider(panel, m_viewer, m_simulator.getNumVerts(), &m_simulator);

		// Re-initialize system and update positions once the user lets go of the slider
		slider->setFinalCallback([this, slider](float v) {
			slider->setValue(v);
			bool wasRunning = m_simActive;
			stop();
			m_simulator.initializeSystem();
			update();
			if (wasRunning) start();
			});

		// For each constraint group add controls over temperature dependent properties
		for (const auto& g : m_simulator.getConstraintGroups()) {

			// Check if the constraint group is temperature dependent (just check the first constraint, they are all of the same type)
			if((g->constraints)[0]->isTemperatureDependent()){

				new Label(m_temperature_window, g->name, "sans-bold");
				Widget* panel = new Widget(m_temperature_window);
				panel->setLayout(new BoxLayout(nanogui::Orientation::Horizontal, nanogui::Alignment::Middle, 0, 10));
				ConstraintPropertySlider* slider = new ConstraintPropertySlider(panel, m_viewer, m_simulator.getNumVerts(), g);

				// Re-initialize system and update positions once the user lets go of the slider
				slider->setFinalCallback([this, slider](float v) {
					slider->setValue(v);
					bool wasRunning = m_simActive;
					stop();
					m_simulator.initializeSystem();
					update();
					if (wasRunning) start();
					});

			}
		}

		m_viewer->performLayout();
	}

    // Called when a new mesh is set in the viewer.
    // We convert the Surface_mesh vertices and triangles into Eigen
    // matrices and pass them to the simulator.
    // Additionally, if desired, we generate tetrahedrons that fill the
    // mesh.
    bool setMesh(bool add_tets) {
        std::cout << "New mesh was loaded, re-initializing simulation..." << std::endl;

        // Stop the running simulation
        stop();

        // Convert surface mesh to Eigen matrices
        surface_mesh::Surface_mesh* mesh = m_viewer->getMesh();
        int j = 0;
        ProjDyn::Positions vertices(mesh->n_vertices(), 3);
        ProjDyn::Triangles faces(mesh->n_faces(), 3);
        for (auto f : mesh->faces()) {
            int k = 0;
            for (auto v : mesh->vertices(f)) {
                faces(j, k) = (ProjDyn::Index)v.idx();
                ++k;
            }
            ++j;
        }
        j = 0;
        for (auto v : mesh->vertices()) {
            vertices.row(j) << (ProjDyn::Scalar)mesh->position(v).x,
                (ProjDyn::Scalar)mesh->position(v).y,
                (ProjDyn::Scalar)mesh->position(v).z;
            ++j;
        }

        ProjDyn::Tetrahedrons tets(0, 4);
        if (add_tets) {
            ProjDyn::Positions vol_verts(0, 3);
            try {
                ProjDyn::tetrahedralizeMesh(vertices, faces, vol_verts, tets, 2);
            }
            catch (int e) {
                std::cout << "Error while generating tet-mesh; error-code: " << e << ". Read the console log above for details." << std::endl;
                return false;
            }
            vertices = vol_verts;
        }

        // Set the mesh in the simulator
        m_simulator.setMesh(vertices, faces, tets);

        // Compute neighbourhood info
        m_vertexStars = ProjDyn::makeVertexStars(vertices.size(), faces);

        if (m_startButton && m_stopButton) {
            m_startButton->setPushed(false);
            m_stopButton->setPushed(true);
        }

        updateConstraintsGUI();
		updateTemperatureGUI();
        m_viewer->showFloor(false);

        return true;
    }

    // Performs a time step and updates the positions that are drawn in the shader window
    bool update(bool forcedUpload = false) {
        if (!m_simulator.isInitialized()) {
            if (!m_simulator.initializeSystem())
                return false;
        }

        // Simulate one time step
        m_simulator.step(m_numIterations);

        return uploadPositions(forcedUpload);
    }

    // Starts a thread that runs the simulation by constantly
    // calling update(), pausing to not run faster than the set FPS
    bool start() {
        stop();

        // Make sure the simulator is properly initialized
        if (!m_simulator.isInitialized()) {
            if (!m_simulator.initializeSystem()) {
                if (m_startButton && m_stopButton) {
                    m_startButton->setPushed(false);
                    m_stopButton->setPushed(true);
                }
                return false;
            }
        }

        // Create a thread that runs the simulation
        // It calls a function that triggers a time-step every 1000/FPS milliseconds
        m_simActive = true;
        m_simulationThread = std::thread(
            [this]() {
            std::chrono::milliseconds time(1000 / FPS);
            while (m_simActive) {
                std::this_thread::sleep_for(time);
                update();
                glfwPostEmptyEvent();
            }
        }
        );

        if (m_startButton && m_stopButton) {
            m_startButton->setPushed(true);
            m_stopButton->setPushed(false);
        }

        return true;
    }

    // Pauses/stops the current simulation by killing the active thread.
    void stop() {
        if (m_simActive) {
            m_simActive = false;
            m_simulationThread.join();
        }
        if (m_startButton && m_stopButton) {
            m_startButton->setPushed(false);
            m_stopButton->setPushed(true);
        }
    }

    // Extract positions, convert them to column-wise tripples of floats and
    // upload them to the OpenGL buffer
    bool uploadPositions(bool forcedUpload = false) {
        const ProjDyn::Positions& pos = m_simulator.getPositions();

        // Initialize matrix if not done already
        if (m_uploadPos.cols() != m_simulator.getNumOuterVerts() || m_uploadPos.rows() != 3) {
            m_uploadPos.resize(3, m_simulator.getNumOuterVerts());
        }

        // There might be more vertices in the simulation, since inner vertices are not present in the m_viewer
#pragma omp parallel for
        for (int i = 0; i < m_uploadPos.cols(); i++) {
            m_uploadPos(0, i) = (float)pos(i, 0);
            m_uploadPos(1, i) = (float)pos(i, 1);
            m_uploadPos(2, i) = (float)pos(i, 2);
        }
        m_viewer->updateShaderVertices(m_uploadPos, forcedUpload);

        if (m_updateNormals && m_vertexStars.size() >= m_simulator.getNumOuterVerts()) {
            // Initialize matrix if not done already
            if (m_uploadNormals.cols() != pos.rows() || m_uploadNormals.rows() != 3) {
                m_uploadNormals.resize(3, m_simulator.getNumOuterVerts());
            }
            // Compute per-triangle normals and sum them on vertices
            const ProjDyn::Triangles& tris = m_simulator.getTriangles();
            m_uploadNormals.setZero();
#pragma omp parallel for
            for (int vInd = 0; vInd < m_simulator.getNumOuterVerts(); vInd++) {
                float fac = 1.f / m_vertexStars[vInd].size();
                ProjDyn::Vector3 normal;
                normal.setZero();
                for (int locInd = 0; locInd < m_vertexStars[vInd].size(); locInd++) {
                    int t = m_vertexStars[vInd][locInd].t1;
                    normal += fac * (pos.row(tris(t, 2)) - pos.row(tris(t, 1))).cross((pos.row(tris(t, 0)) - pos.row(tris(t, 1)))).normalized();
                }
                m_uploadNormals(0, vInd) = normal(0);
                m_uploadNormals(1, vInd) = normal(1);
                m_uploadNormals(2, vInd) = normal(2);
            }

			// Upload the colors
			updateColor(); //New (to change color)

            //Upload the normals
            m_viewer->updateShaderNormals(m_uploadNormals, forcedUpload);
        }

        return true;
    }

	void updateColor() {
		ProjDyn::Vector temperatures = m_simulator.getTemperatures();
		ProjDyn::Index n_vertices = m_simulator.getNumVerts();
		MatrixXf vertex_color(3, n_vertices);
		vertex_color.setZero();
		for (int i = 0; i < n_vertices; i++) vertex_color.col(i) = Vector3f(temperatures[i]/100, 0.2, 1.0 - temperatures[i] / 100);
		m_viewer->changeColor(vertex_color);
		std::cout << "Mean temperature: " << temperatures.mean() << std::endl;
	}

    // Set external forces to point into downwards y direction with a certain magnitude
    void setGravity(ProjDyn::Scalar g) {
        m_simulator.setGravity(g);
    }

    // Can be used to visually mark certain vertices based on an integer
    // value they receive in a row-vector (vertex ID corresponds to column)
    void uploadVertexStatus(const Eigen::Matrix<int, 1, -1>& vStatus) {
        m_viewer->updateVertexStatus(vStatus);
    }

    // Change the number of local-global iterations in the simulator
    void setNumIterations(Index numIts) {
        m_numIterations = numIts;
    }
    Index getNumIterations() {
        return m_numIterations;
    }

    // Called when a new mesh was set in the viewer
    bool setMesh() {
        return setMesh(false);
    }

    // Add constraints to the simulator, either as groups, as lists
    // or as single constraints
    void addConstraints(const std::vector<ProjDyn::ConstraintPtr>& constraints) {
        m_simulator.addConstraints(constraints);
        updateConstraintsGUI();
		updateTemperatureGUI();
    }
    void addConstraints(ProjDyn::ConstraintGroupPtr constraints) {
        m_simulator.addConstraints(constraints);
        updateConstraintsGUI();
		updateTemperatureGUI();
    }
    void addConstraint(const ProjDyn::ConstraintPtr& constraint) {
        m_simulator.addConstraint(constraint);
        updateConstraintsGUI();
		updateTemperatureGUI();
    }

    // Remove a constraint group from the simulator, either by value or by name.
    void removeConstraints(ProjDyn::ConstraintGroupPtr constraints) {
        m_simulator.removeConstraints(constraints);
        updateConstraintsGUI();
		updateTemperatureGUI();
    }
    void removeConstraints(std::string constraint_name) {
        m_simulator.removeConstraints(constraint_name);
        updateConstraintsGUI();
		updateTemperatureGUI();
    }

    // Remove all constraints from the simulation.
    void clearConstraints() {
        m_simulator.clearConstraints();
    }

    // In the following are several helper function to add several types of specific
    // constraints to the simulation:


    // Add tetrahedral strain constraints to all tets:
    void addTetStrainConstraints(ProjDyn::Scalar weight = 1.) {
        std::vector<Index> allTets;
        const ProjDyn::Tetrahedrons& tets = m_simulator.getTetrahedrons();
        for (Index i = 0; i < tets.rows(); i++) allTets.push_back(i);
        addTetStrainConstraints(allTets, weight);
    }

    // Add tetrahedral strain constraints to some tets:
    void addTetStrainConstraints(const std::vector<Index>& tetInds, ProjDyn::Scalar weight = 1.) {
        std::vector<ProjDyn::ConstraintPtr> tet_constraints;
        const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
        const ProjDyn::Tetrahedrons& tets = m_simulator.getTetrahedrons();
        if (tets.rows() > 0) {
            for (Index i : tetInds) {
                if (i >= tets.rows()) continue;
                std::vector<ProjDyn::Index> tetInds;
                for (int j = 0; j < 4; j++) tetInds.push_back(tets(i, j));
                // The weight is the tet volume
                ProjDyn::Scalar w = ProjDyn::tetrahedronVolume(sim_verts, tets.row(i));
                if (w > 1e-6) {
                    // The constraint is constructed, made into a shared pointer and appended to the list
                    // of constraints.
                    ProjDyn::TetStrainConstraint* esc = new ProjDyn::TetStrainConstraint(tetInds, w, sim_verts);
                    tet_constraints.push_back(std::shared_ptr<ProjDyn::TetStrainConstraint>(esc));
                }
            }
            addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Tet Strain", tet_constraints, weight));
        }
        else {
            std::cout << "WARNING: Could not add tet-strain - no tets available!" << std::endl;
        }
    }

    // Add bending constraints to all vertices of the simulation
    void addBendingConstraints(ProjDyn::Scalar weight = 1.) {
        std::vector<Index> allOuterVerts;
        for (Index i = 0; i < m_simulator.getNumOuterVerts(); i++) allOuterVerts.push_back(i);
        addBendingConstraints(allOuterVerts, weight);
    }

    // Add bending constraints to some vertices of the simulation
    void addBendingConstraints(const std::vector<Index>& vertInds, ProjDyn::Scalar weight = 1.) {
        std::vector<ProjDyn::ConstraintPtr> bend_constraints;
        const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
        std::vector<ProjDyn::VertexStar>& vStars = m_vertexStars;
        ProjDyn::Vector voronoiAreas = ProjDyn::vertexMasses(m_simulator.getPositions(), m_simulator.getTriangles());
        Surface_mesh* smesh = m_viewer->getMesh();
        for (auto v : smesh->vertices()) {
            Index i = v.idx();
            if (std::find(vertInds.begin(), vertInds.end(), i) == vertInds.end()) continue;
            if (smesh->is_boundary(v)) continue;
            if (i >= m_simulator.getNumOuterVerts()) continue;
            // The weight is the voronoi area
            ProjDyn::Scalar w = voronoiAreas(i) * 0.01;
            if (w > 1e-6) {
                // The constraint is constructed, made into a shared pointer and appended to the list
                // of constraints.
                ProjDyn::BendingConstraint* esc = new ProjDyn::BendingConstraint(vStars[i], w, voronoiAreas(i), sim_verts, m_simulator.getTriangles());
                bend_constraints.push_back(std::shared_ptr<ProjDyn::BendingConstraint>(esc));
            }
        }
        addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Tri Bending", bend_constraints, weight));
    }

    // Add triangular strain constraints to all triangles of the simulation
    void addTriangleStrainConstraints(ProjDyn::Scalar weight = 1.) {
        std::vector<Index> allTris;
        for (Index i = 0; i < m_simulator.getTriangles().rows(); i++) allTris.push_back(i);
        addTriangleStrainConstraints(allTris, weight);
    }

    // Add triangular strain constraints to some triangles of the simulation
    void addTriangleStrainConstraints(const std::vector<Index>& triInds, ProjDyn::Scalar weight = 1.) {
        std::vector<ProjDyn::ConstraintPtr> tri_constraints;
        const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
        const ProjDyn::Triangles& tris = m_simulator.getTriangles();
        for (Index i : triInds) {
            if (i >= tris.rows()) continue;
            std::vector<ProjDyn::Index> triInds;
            for (int j = 0; j < 3; j++) triInds.push_back(tris(i, j));
            // The weight is the triangle area
            ProjDyn::Scalar w = ProjDyn::triangleArea(sim_verts, tris.row(i));
            if (w > 1e-6) {
                // The constraint is constructed, made into a shared pointer and appended to the list
                // of constraints.
                ProjDyn::TriangleStrainConstraint* esc = new ProjDyn::TriangleStrainConstraint(triInds, w, sim_verts);
                tri_constraints.push_back(std::shared_ptr<ProjDyn::TriangleStrainConstraint>(esc));
            }
        }
        // Add constraints
        addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Tri Strain", tri_constraints, weight));
    }

    // Add edge spring constraints to all edges of the simulation
    void addEdgeSpringConstraints(ProjDyn::Scalar weight = 1.) {
        // For tet meshes we cannot use Surface_mesh
        if (m_simulator.getTetrahedrons().rows() > 0) {
            addEdgeSpringConstraintsTets(weight);
        }
        else {
            const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
            const ProjDyn::Triangles& tris = m_simulator.getTriangles();
            Surface_mesh* smesh = m_viewer->getMesh();
            std::vector<ProjDyn::ConstraintPtr> spring_constraints;
            for (auto edge : smesh->edges()) {
                // The weight is set to the edge length
                ProjDyn::Scalar w = (sim_verts.row(smesh->vertex(edge, 0).idx()) - sim_verts.row(smesh->vertex(edge, 1).idx())).norm();
                if (w > 1e-6) {
                    // The constraint is constructed, made into a shared pointer and appended to the list
                    // of constraints.
                    std::vector<Index> edge_inds;
                    edge_inds.push_back(smesh->vertex(edge, 0).idx());
                    edge_inds.push_back(smesh->vertex(edge, 1).idx());
                    ProjDyn::EdgeSpringConstraint* esc = new ProjDyn::EdgeSpringConstraint(edge_inds, w, sim_verts);
                    spring_constraints.push_back(std::shared_ptr<ProjDyn::EdgeSpringConstraint>(esc));
                }
            }
            addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Edge Springs", spring_constraints, weight));
        }
    }

    // Gets called when the user is grabbing some vertices with the mouse
    // (alt + left-click).
    void grab(const std::vector<ProjDyn::Index>& grabbedVerts, const std::vector<Vector3f>& grabPos) {
        if (!m_simulator.isInitialized() || !m_simActive) return;

        m_simulator.setGrab(grabbedVerts, grabPos);
    }

    // Gets called when the user is letting go of the grabbed vertices.
    void releaseGrab() {
        if (!m_simulator.isInitialized()) return;

        m_simulator.releaseGrab();
    }

    // In the following are some simple getter/setter function to communicate
    // with the simulation:

    const ProjDyn::Positions& getPositions() {
        return m_simulator.getPositions();
    }

    const ProjDyn::Triangles& getTriangles() {
        return m_simulator.getTriangles();
    }

    const ProjDyn::Tetrahedron& getTetrahedra() {
        return m_simulator.getTetrahedrons();
    }

    const ProjDyn::Index getNumVertices() {
        return m_simulator.getNumVerts();
    }

private:
    // Global variables used by the callback functions
    int m_numIterations = NUM_ITS_INITIAL;
    ProjDyn::Simulator m_simulator;
    Eigen::MatrixXf m_uploadPos, m_uploadNormals;
    std::thread m_simulationThread;
    bool m_simActive = false;
    std::vector<ProjDyn::VertexStar> m_vertexStars;
    Window* m_constraint_window = nullptr;
	Window* m_temperature_window = nullptr;
    Viewer* m_viewer;
    Button* m_startButton = nullptr;
    Button* m_stopButton = nullptr;
    bool m_updateNormals = UPDATE_NORMALS;

    void addEdgeSpringConstraintsTets(ProjDyn::Scalar weight = 1.) {
        const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
        std::vector<ProjDyn::ConstraintPtr> spring_constraints;
        const ProjDyn::Tetrahedrons& tets = m_simulator.getTetrahedrons();
        // If tets are available, add a spring on each tet-edge
        for (Index i = 0; i < tets.rows(); i++) {
            for (int j = 0; j < 4; j++) {
                std::vector<ProjDyn::Index> edge;
                edge.push_back(tets(i, j));
                edge.push_back(tets(i, (j + 1) % 4));
                // Easy way to make sure each edge only gets added once:
                if (edge[0] < edge[1]) {
                    // The weight is set to the edge length
                    ProjDyn::Scalar w = (sim_verts.row(edge[0]) - sim_verts.row(edge[1])).norm();
                    if (w > 1e-6) {
                        // The constraint is constructed, made into a shared pointer and appended to the list
                        // of constraints.
                        ProjDyn::EdgeSpringConstraint* esc = new ProjDyn::EdgeSpringConstraint(edge, w, sim_verts);
                        spring_constraints.push_back(std::shared_ptr<ProjDyn::EdgeSpringConstraint>(esc));
                    }
                }
            }
        }

        // Add constraints
        addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Edge Springs", spring_constraints, weight));
    }
};