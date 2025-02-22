//
// Created by jan on 11/29/22.
//

#include "io/input/sim_input/XMLReader.h"


namespace io::input {
    XMLReader::XMLReader() = default;

    XMLReader::~XMLReader() = default;

    void XMLReader::readFile(const char *filename, std::list<Particle> &particles, std::list<Membrane>& membranes,
                             std::unordered_map<io::input::names, std::string> &arg_map) {
        xml_schema::properties properties;
        properties.no_namespace_schema_location("./XMLFormat.xsd");

        try {
            std::unique_ptr<simulation_t> simulation {Simulation(std::string{filename}, 0, properties)};
            // we will prefer arguments from a checkpoint file to a normal one
            bool isCheckpoint = simulation->FileType().Checkpoint().present();
            /**Set str in arg map at name if present and if arg map has no value there or if this file is a checkpoint file*/
            auto setInMap = [&](names name, bool present, const std::string& def, const std::function<std::string()>& str_get){
                if (!present && !arg_map.contains(name)) {
                    arg_map.emplace(name, def);
                    return;
                }
                if (isCheckpoint || !arg_map.contains(name)) arg_map.emplace(name, str_get());
            };
            /**Set str in arg map at name if arg map has no value there or if this file is a checkpoint file*/
            auto setInMapND = [&](names name, const std::string& str) {
                if (isCheckpoint || !arg_map.contains(name)) arg_map.emplace(name, str);
            };
            std::string parseBuffer;
            auto setInMapVal = [&](names name, bool present, const std::string& def, const std::function<std::string()>& str_get){
                if (!present && !arg_map.contains(name)) {
                    arg_map.emplace(name, def);
                    parseBuffer = def;
                    return;
                }
                if (isCheckpoint || !arg_map.contains(name)) {
                    arg_map.emplace(name, str_get());
                    parseBuffer = str_get();
                }
                parseBuffer = arg_map.at(name);
            };


            // <!-- IO -->

            if(simulation->OutputFile().present()) {
                setInMap(outputFilePath, simulation->OutputFile()->FolderPath().present(), default_output_folder, [&]()->std::string{return simulation->OutputFile()->FolderPath().get();});
                setInMap(outputFileName, simulation->OutputFile()->OutputFileName().present(), default_output_base_name, [&]()->std::string{return simulation->OutputFile()->OutputFileName().get();});
            }

            // <!-- Calculation Decisions -->

            setInMap(startTime, simulation->StartTime().present(), std::to_string(default_start_time), [&]()->std::string{return std::to_string(simulation->StartTime().get());});
            setInMap(endTime, simulation->EndTime().present(), std::to_string(default_end_time), [&]()->std::string{return std::to_string(simulation->EndTime().get());});
            setInMap(delta_t, simulation->TimeStepSize().present(), std::to_string(default_delta_t), [&]()->std::string{return std::to_string(simulation->TimeStepSize().get());});

            if (simulation->ForceCalculation().Gravity().present()) {
                setInMapND(forceCalculation, "gravity");
            }
            else if (auto& lj = simulation->ForceCalculation().LennardJones(); lj.present()) {
                setInMapND(forceCalculation, "lennardJones");
                setInMap(epsilon, lj->Epsilon().present(), std::to_string(default_epsilon), [&]()->std::string{return std::to_string(lj->Epsilon().get());});
                setInMap(sigma, lj->Sigma().present(), std::to_string(default_sigma), [&]()->std::string{return std::to_string(lj->Sigma().get());});
            }
            else {
                output::loggers::general->debug("This really shouldn't happen. No ForceCalculation was specified despite it being mandatory. Using default...");
                setInMapND(forceCalculation, default_force_type);
            }

            if (auto& eGrav = simulation->ForceCalculation().EnableGrav(); eGrav.present()) {
                setInMapND(enableGrav, "1");
                setInMapND(gGrav0, std::to_string(eGrav->X()));
                setInMapND(gGrav1, std::to_string(eGrav->Y()));
                setInMapND(gGrav2, std::to_string(eGrav->Z()));
            } else {
                setInMapND(enableGrav, "0");
            }
            if (auto& eLC = simulation->ForceCalculation().EnableLC(); eLC.present()) {
                setInMapND(enableLinkedCell, "1");
                setInMapND(rCutoff, std::to_string(eLC->CutoffRadius()));
                setInMapND(boundingBox_X0, std::to_string(eLC->BoundaryBox().BoxSize().X()));
                setInMapND(boundingBox_X1, std::to_string(eLC->BoundaryBox().BoxSize().Y()));
                setInMapND(boundingBox_X2, std::to_string(eLC->BoundaryBox().BoxSize().Z()));

                setInMapND(boundCondFront, eLC->BoundaryBox().Front());
                setInMapND(boundCondRear, eLC->BoundaryBox().Rear());
                setInMapND(boundCondLeft, eLC->BoundaryBox().Left());
                setInMapND(boundCondRight, eLC->BoundaryBox().Right());
                setInMapND(boundCondTop, eLC->BoundaryBox().Top());
                setInMapND(boundCondBottom, eLC->BoundaryBox().Bottom());

                setInMap(enableProfiler, eLC->EnableProfiling().present(), std::to_string(0), [&]()->std::string{return std::to_string(1);});
                if(eLC->EnableProfiling().present()){
                    setInMapND(profilerNumBins, std::to_string(eLC->EnableProfiling()->NumBins()));
                }
            } else {
                setInMapND(enableLinkedCell, "0");
            }
            setInMapND(enableOMP, std::to_string(static_cast<int>(simulation->ForceCalculation().EnableOMP().present())));
            if (auto& eMem = simulation->ForceCalculation().EnableMem(); eMem.present()) {
                setInMapND(enableMembrane, "1");
                setInMapND(enableMembranePull, std::to_string(static_cast<int>(eMem->EnableMemPull().present())));
            } else {
                setInMapND(enableMembrane, "0");
                setInMapND(enableMembranePull, "0");
            }

            setInMap(positionCalculation, simulation->PositionCalculation().present(), default_pos_type,[&]()->std::string{return simulation->PositionCalculation().get();});
            setInMap(velocityCalculation, simulation->VelocityCalculation().present(), default_vel_type,[&]()->std::string{return simulation->VelocityCalculation().get();});

            setInMapVal(brown, simulation->AverageBrownianMotion().present(), std::to_string(default_brown), [&]()->std::string{return std::to_string(simulation->AverageBrownianMotion().get());});
            double brown_val = std::stod(parseBuffer);

            setInMapVal(dimensions, simulation->Dimensions().present(), std::to_string(default_dims), [&]()->std::string{return std::to_string(simulation->Dimensions().get());});
            int dims_val = std::stoi(parseBuffer);

            if (auto& t = simulation->Thermostat(); t.present()) {
                setInMapND(enableThermo, std::to_string(1));
                setInMapND(thermoTInit, std::to_string(t.get().T_Init()));
                setInMapND(thermoNTerm, std::to_string(t.get().N_Term()));
                setInMap(thermoTTarget, t.get().T_Target().present(), std::to_string(t.get().T_Init()), [&]()->std::string{return std::to_string(t.get().T_Target().get());});
                setInMap(thermoDelta_t, t.get().Delta_T().present(), std::to_string(default_delta_temp), [&]()->std::string{return std::to_string(t.get().Delta_T().get());});
                setInMap(thermoType_t, t.get().ThermoMode().present(), std::to_string(default_thermo_mode), [&]()->std::string{return std::to_string(t.get().ThermoMode().get() == "Normal"?0:1);});
            } else {
                setInMapND(enableThermo, std::to_string(0));
            }


            // <!-- Misc -->

            setInMap(logLevel, simulation->LogLevel().present(), std::to_string(default_log_level), [&]()->std::string{return std::to_string(simulation->LogLevel().get());});
            setInMap(enableCheckpointing, simulation->EnableCheckpointing().present(), std::to_string(default_checkpointing), [&]()->std::string{return std::to_string(simulation->EnableCheckpointing().get());});

            if (simulation->Benchmark().present()) {
                setInMapND(benchmark, "1");
                if (simulation->Benchmark()->BenchmarkType().FileBenchmark().present()) {
                    setInMapND(benchmarkType, "file");
                }
                else {
                    setInMapND(benchmarkType, "default");
                    if (simulation->Benchmark()->BenchmarkType().DefaultBenchmark()->MaximumBodySize().present()) {
                        setInMapND(benchMaxBodySize, std::to_string(simulation->Benchmark()->BenchmarkType().DefaultBenchmark()->MaximumBodySize().get()));
                    }
                    else {
                        setInMapND(benchMaxBodySize, std::to_string(default_bench_maxBody));
                    }
                }
                setInMap(benchIterationCount, simulation->Benchmark()->IterationCount().present(), std::to_string(default_bench_iterations), [&]()->std::string{return std::to_string(simulation->Benchmark()->IterationCount().get());});
            }

            //handle body/particles according to file type
            if (simulation->FileType().Input().present()) {
                auto& shapeList = simulation->FileType().Input().get().ShapeList();
                for (const auto & s : shapeList.Shape()) {
                    Body body;

                    if (s.Particle().present()) {
                        body.shape = Shape::particle;

                        dvectorToEigenVector3d(s.Particle()->Position(), body.fixpoint);
                        dvectorToEigenVector3d(s.Particle()->Velocity(), body.start_velocity);

                        if (s.Particle()->Mass() == 0) {
                            output::loggers::general->warn("Particle has a mass of 0, which is illegal. Skipping this particle...");
                            continue;
                        }

                        body.mass = s.Particle()->Mass();
                        double eps, sig;
                        if (s.Particle()->Sigma().present()) sig = s.Particle()->Sigma().get();
                        else if (arg_map.contains(sigma)) sig = std::stod(arg_map[sigma]);
                        else sig = default_sigma;
                        if (s.Particle()->Epsilon().present()) eps = s.Particle()->Epsilon().get();
                        else if (arg_map.contains(epsilon)) eps = std::stod(arg_map[epsilon]);
                        else eps = default_epsilon;
                        ParticleGenerator::generateParticle(body.fixpoint, body.start_velocity, body.mass, particles, sig, eps);
                    }
                    else if (s.Cuboid().present()) {
                        body.shape = Shape::cuboid;

                        dvectorToEigenVector3d(s.Cuboid()->Position(), body.fixpoint);
                        dvectorToEigenVector3d(s.Cuboid()->Velocity(), body.start_velocity);

                        ivectorToEigenVector3d(s.Cuboid()->Dimensions(), body.dimensions);

                        if (s.Cuboid()->Mass() == 0 || s.Cuboid()->Spacing() == 0) {
                            output::loggers::general->warn("Cuboid has a mass or spacing of 0, which is illegal. Skipping this cuboid...");
                            continue;
                        }
                        body.distance = s.Cuboid()->Spacing();
                        body.mass = s.Cuboid()->Mass();

                        double eps, sig;
                        if (s.Cuboid()->Sigma().present()) sig = s.Cuboid()->Sigma().get();
                        else if (arg_map.contains(sigma)) sig = std::stod(arg_map[sigma]);
                        else sig = default_sigma;
                        if (s.Cuboid()->Epsilon().present()) eps = s.Cuboid()->Epsilon().get();
                        else if (arg_map.contains(epsilon)) eps = std::stod(arg_map[epsilon]);
                        else eps = default_epsilon;
                        ParticleGenerator::generateCuboid(body, brown_val, particles, dims_val, sig, eps);
                    }

                    else if (s.Sphere().present()) {
                        body.shape = Shape::sphere;

                        dvectorToEigenVector3d(s.Sphere()->Position(), body.fixpoint);
                        dvectorToEigenVector3d(s.Sphere()->Velocity(), body.start_velocity);

                        body.dimensions << s.Sphere()->Radius(), s.Sphere()->Radius(), s.Sphere()->Radius();

                        if (s.Sphere()->Mass() == 0 || s.Sphere()->Spacing() == 0) {
                            output::loggers::general->warn("Sphere has a mass or spacing of 0, which is illegal. Skipping this sphere...");
                            continue;
                        }
                        body.distance = s.Sphere()->Spacing();
                        body.mass = s.Sphere()->Mass();
                        double eps, sig;
                        if (s.Sphere()->Sigma().present()) sig = s.Sphere()->Sigma().get();
                        else if (arg_map.contains(sigma)) sig = std::stod(arg_map[sigma]);
                        else sig = default_sigma;
                        if (s.Sphere()->Epsilon().present()) eps = s.Sphere()->Epsilon().get();
                        else if (arg_map.contains(epsilon)) eps = std::stod(arg_map[epsilon]);
                        else eps = default_epsilon;
                        ParticleGenerator::generateSphere(body, brown_val, particles, dims_val, sig, eps);
                    }

                    //How do i add the membrane as an option here? ShapeList is automatically generated..
                    else if (s.Membrane().present()) {
                        body.shape = Shape::membrane;

                        dvectorToEigenVector3d(s.Membrane()->Position(), body.fixpoint);
                        dvectorToEigenVector3d(s.Membrane()->Velocity(), body.start_velocity);
                        ivectorToEigenVector3d(s.Membrane()->Dimensions(), body.dimensions);

                        if (auto& pull = s.Membrane()->Pull(); pull.present()) {
                            body.pullEndTime = (pull->EndTimePull().present()) ? pull->EndTimePull().get() : std::stod(arg_map[endTime]);
                            body.pullForce = { pull->PullForce().X(), pull->PullForce().Y(), pull->PullForce().Z() };
                            for(auto& ind : pull->PullIndices().Index()) body.pullIndices.push_back({ind.I(), ind.J()});
                        } else {
                            body.pullEndTime = std::stod(arg_map[startTime]);
                            body.pullForce = { 0 ,0 ,0 };
                            body.pullIndices.clear();
                        }

                        if (s.Membrane()->Mass() == 0 || s.Membrane()->Spacing() == 0) {
                            output::loggers::general->warn("Membrane has a mass or spacing of 0, which is illegal. Skipping this membrane...");
                            continue;
                        }
                        body.distance = s.Membrane()->Spacing();
                        body.mass = s.Membrane()->Mass();
                        double eps, sig;
                        if (s.Membrane()->Sigma().present()) sig = s.Membrane()->Sigma().get();
                        else if (arg_map.contains(sigma)) sig = std::stod(arg_map[sigma]);
                        else sig = default_sigma;
                        if (s.Membrane()->Epsilon().present()) eps = s.Membrane()->Epsilon().get();
                        else if (arg_map.contains(epsilon)) eps = std::stod(arg_map[epsilon]);
                        else eps = default_epsilon;

                        body.desiredDistance = s.Membrane()->DesiredDistance();
                        body.springStrength = s.Membrane()->SpringStrength();
                        ParticleGenerator::generateMembrane(body, brown_val, particles, membranes, dims_val, sig, eps);
                    }
                    else {
                        output::loggers::general->debug("An unknown shape was detected. This really shouldn't happen. Skipping...");
                        continue;
                    }

                    // handle body and create particles

                }
            }
            else if (simulation->FileType().Checkpoint().present()) {
                setInMapND(simLastIteration, std::to_string(simulation->FileType().Checkpoint().get().LastIteration()));

                auto& cpParticles = simulation->FileType().Checkpoint().get().CPParticles();
                for (auto& cp : cpParticles.CPParticle()) {
                    Particle p {std::array<double, 3>{cp.Position().X(), cp.Position().Y(), cp.Position().Z()},
                                std::array<double, 3>{cp.Velocity().X(), cp.Velocity().Y(), cp.Velocity().Z()},
                                cp.Mass(), cp.Type()};
                    p.setF({cp.Force().X(), cp.Force().Y(), cp.Force().Z()});
                    p.setOldF({cp.OldForce().X(), cp.OldForce().Y(), cp.OldForce().Z()});
                    p.setSigma(cp.Sigma());
                    p.setEpsilon(cp.Epsilon());
                    p.setID(cp.Id());
                    particles.push_back(p);
                }
            }
            else {
                output::loggers::general->error("unsupported FileType in XML. Valid are Input or Checkpoint");
                exit(-1);
            }


        }
        catch (const xml_schema::exception &e) {
            output::loggers::general->error("The following exception occurred during the parsing of your XML input file:");
            output::loggers::general->error(e.what());
            exit(1);
        }
    }

    void XMLReader::dvectorToEigenVector3d(dvector_t const &dv, Eigen::Vector3d &ev) {
        ev << dv.X(), dv.Y(), dv.Z();
    }

    void XMLReader::ivectorToEigenVector3d(ivector_t const &dv, Eigen::Vector3d &ev) {
        ev << dv.X(), dv.Y(), dv.Z();
    }
}