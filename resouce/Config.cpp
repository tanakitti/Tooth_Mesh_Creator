//#include "Config.h"
//
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <iterator>
//#include <iostream>
//
//namespace pd
//{
//
//    Config::Config()
//    {
//        std::cout << "[Config] Reading Config file: " << ConfigFileName << std::endl;
//
//        readParameterBool(  &m_config.collisionsLimitFps,                                    "collisionsLimitFps",                                     1, "Whether to limit the collision thread's execution frequency." );
//        readParameterInt(   &m_config.collisionsFpsTarget,                                   "collisionsFpsTarget",                                    1, "The frequency that the collision thread's execution should not exceed." );
//        readParameterInt(   &m_config.collisionsFpsUpdateFreq,                               "collisionsFpsUpdateFreq",                                1, "The frequency (per second) with which the current framerate of the collision thread is estimated." );
//
//        readParameterBool(  &m_config.simulationUseMultipleGpus,                             "simulationUseMultipleGpus",                              1, "Whether to use two CUDA devices. If yes, the following device Id's will be switched at runtime.\n# Otherwise, the default CUDA device is used for everything." );
//        readParameterInt(   &m_config.simulationComputeGpuId,                                "simulationComputeGpuId",                                 1, "CUDA id of the primary device for simulation." );
//        readParameterInt(   &m_config.simulationVisualizeGpuId,                              "simulationVisualizeGpuId",                               1, "CUDA id of the secondary device for simulation, will do scalar field generation." );
//        readParameterFloat(  m_config.simulationInDeviceTorqueMods,                          "simulationInDeviceTorqueMods",                           3, "CUDA id of the secondary device for simulation, will do scalar field generation." );
//
//        readParameterFloat( &m_config.simulationFrictionNormalForceMod,                      "simulationFrictionNormalForceMod",                       1, "Friction material coefficient (coulomb model) when the tool is in sliding contact. (Higher value => more friction when sliding)" );
//        readParameterFloat( &m_config.simulationFrictionKinecticNormalForceSlidingThreshold, "simulationFrictionKinecticNormalForceSlidingThreshold",  1, "Friction material coefficient (coulomb model) when the tool is in sliding contact. (Higher value => more friction when sliding)" );
//        readParameterFloat( &m_config.simulationFrictionStaticNormalForceSlidingThreshold,   "simulationFrictionStaticNormalForceSlidingThreshold",    1, "Friction material coefficient (coulomb model) when the tool is in static contact. (Higher value => more friction when in static contact)" );
//        readParameterFloat( &m_config.simulationFrictionSlidingPushback,                     "simulationFrictionSlidingPushback",                      1, "The distance the tool is pushed back from the surface to move laterally. (High value => less friction) Should only be adjusted manually once, so that the friction is in the correct neighborhood." );
//        readParameterFloat( &m_config.simulationFrictionDrillingSlidingMod,                  "simulationFrictionDrillingSlidingMod",                   1, "A modifier that is designed to reduce friction when the tool is drilling. (Smaller value => more friction when drilling)" );
//        readParameterFloat( &m_config.simulationContactEpsilon,                              "simulationContactEpsilon",                               1, "A small margin of error that is acceptable to travel backwards after a contact." );
//        readParameterFloat( &m_config.simulationContactRadiusMod,                            "simulationContactRadiusMod",                             1, "Ratio by which the radius of the tool sphere is multiplied to calculate the penalty force at contact with.\n# Needs to be larger than 1 (>!1)" );
//
//        readParameterFloat( &m_config.simulationDrillRadiusMod,                              "simulationDrillRadiusMod",                               1, "Whether sphere's should be enlarged when drilling (high value => larger drill result)" );
//        readParameterFloat( &m_config.simulationDrillSpeedMod,                               "simulationDrillSpeedMod",                                1, "Drill speed per frame (this is a ratio of depth of penetration, so the deeper the penetration, the faster it will drill), should be very small." );
//        readParameterFloat( &m_config.simulationDrillSpeedMin,                               "simulationDrillSpeedMin",                                1, "Depth minimum after which drilling happens." );
//        readParameterFloat( &m_config.simulationDrillSpeedMax,                               "simulationDrillSpeedMax",                                1, "Maximum drill speed that can be achieved by the dop ratio simulationDrillSpeedMod (should be very small for hip reaming)." );
//        readParameterBool(  &m_config.simulationDrillOff,                                    "simulationDrillOff",                                     1, "Whether to never drill (even if the button is pressed)." );
//
//        readParameterFloat( &m_config.simulationDensityContactVolumeMod,                     "simulationDensityContactVolumeMod",                      1, "Density-scaled sum of contact volume is scaled with this before speed reduction. (Higher value => slower drilling with same contact surface)" );
//        readParameterFloat( &m_config.simulationDensityNormalizedContactVolumeMod,           "simulationDensityNormalizedContactVolumeMod",            1, "Normalized density-scaled sum of contact volume is scaled by this. (Higher value => slower drilling with same contact surface) Adjust this before simulationDensityContactVolumeMod." );
//        readParameterFloat( &m_config.simulationDensityDampingMod,                           "simulationDensityDampingMod",                            1, "Decrease the effect differently dense spheres do for the force and normal density-weighted sum do. (Smaller value => density is less important)" );
//
//        readParameterFloat( &m_config.visualizationDefaultScalarfieldValue,                  "visualizationDefaultScalarfieldValue",                   1, "The default value that all scalar field cells are initialized with before filling with sphere surface." );
//        readParameterInt( (int*)&m_config.visualizationSimilarityMetric,                     "visualizationSimilarityMetric",                          1, "Choose which similarity metric to use to compare user's drill result with ideal drill result. (0=Dice, 1=Jaccard, 2=Accuracy, 3=Sensitivity_TPR, 4=Specificity, 5=Fallout_FPR, 6=AUC, 7=Precision, 8=Recall, 9=CohenKappa, 10=FalsePositive)" );
//        readParameterBool(  &m_config.visualizationSmoothScalarfield,                        "visualizationSmoothScalarfield",                         1, "Whether to smooth the isosurface." );
//        readParameterBool(  &m_config.visualizationSmoothScalarfieldOnGpu,                   "visualizationSmoothScalarfieldOnGpu",                    1, "Whether to smooth the isosurface on the GPU. (There currently is no CPU option, so leave it on)" );
//        readParameterInt( (int*)&m_config.visualizationDrillComparisonSmoothMethod,          "visualizationDrillComparisonSmoothMethod",               1, "The smoothing method that is used after calculating the volume difference and intersection between user drill result and ideal drill result. (0=simple, 1=smoothing #1, 2=smoothing #2 (smoothing #2 gives prettiest results in general))" );
//        readParameterInt(   &m_config.visualizationSmoothGridRes,                            "visualizationSmoothGridRes",                             1, "Dimension of the neighbourhood around a scalar value to smooth the surface.\n# (Should be an odd number, otherwise it's assumed to be N+1)" );
//        readParameterFloat( &m_config.visualizationSmoothSigmaSimilarity,                    "visualizationSmoothSigmaSimilarity",                     1, "Controls how similar neighbours have to be to affect each other's smoothing result. (high value => more smoothing)" );
//        readParameterFloat( &m_config.visualizationSmoothSigmaDistance,                      "visualizationSmoothSigmaDistance",                       1, "Controls how close neighbours have to be to affect each other's smoothing result. (high value => more smoothing)" );
//        readParameterBool(  &m_config.visualizationAdjustScore,                              "visualizationAdjustScore",                               1, "Whether to adjust the score-range based on the following parameters." );
//        readParameterFloat( &m_config.visualizationScoreNeutralSimilarity,                   "visualizationScoreNeutralSimilarity",                    1, "The the default similarity measurement, that will be scored with visualizationScoreNeutralScore." );
//        readParameterFloat( &m_config.visualizationScoreNeutralScore,                        "visualizationScoreNeutralScore",                         1, "The desired similarity function's score, when no material was removed." );
//
//        readParameterFloat( &m_config.hapticsForceMod,                                       "hapticsForceMod",                                        1, "Linear spring hardness. Strength of the spring between haptic device position and surface position. (High value => more force for same dop)\n# Does not work with Non-Linear Growing Force." );
//        readParameterFloat( &m_config.hapticsTorqueMod,                                      "hapticsTorqueMod",                                       1, "Angular spring hardness. Calculated torque is multiplied by this. (High value => more torque)" );
//        readParameterBool(  &m_config.hapticsLimitFps,                                       "hapticsLimitFps",                                        1, "Whether to limit haptic thread execution speed." );
//        readParameterInt(   &m_config.hapticsFpsTarget,                                      "hapticsFpsTarget",                                       1, "The frequency that the haptic thread should not exceed." );
//        readParameterInt(   &m_config.hapticsFpsUpdateFreq,                                  "hapticsFpsUpdateFreq",                                   1, "The frequency (per second) with which the current framerate of the haptic thread is estimated." );
//        readParameterFloat( &m_config.hapticsNonKukaDampingLinear,                           "hapticsNonKukaDampingLinear",                            1, "The linear damping (as a ratio of linear velocity) that non-KUKA haptic devices are dampened by." );
//        readParameterInt(   &m_config.hapticsPoseHistorySlots,                               "hapticsPoseHistorySlots",                                1, "The amount of frames to keep as a running history to estimate the current velocity." );
//        readParameterInt(   &m_config.hapticsTimeDeltaLinearVelocityEstimation,              "hapticsTimeDeltaLinearVelocityEstimation",               1, "The time delta (in ms) between this frame and the history that is used to estimate the linear velocity. (high value => less dynamic but slower updating estimation)" );
//        readParameterInt(   &m_config.hapticsTimeDeltaAngularVelocityEstimation,             "hapticsTimeDeltaAngularVelocityEstimation",              1, "The time delta (in ms) between this frame and the history that is used to estimate the angular velocity. (high value => less dynamic but slower updating estimation)" );
//        readParameterFloat( &m_config.hapticsLinearVelocityLimit,                            "hapticsLinearVelocityLimit",                             1, "The linear velocity limit (in m/s) that the haptic device should not exceed." );
//        readParameterFloat( &m_config.hapticsAngularVelocityLimit,                           "hapticsAngularVelocityLimit",                            1, "The angular velocity limit (in degree/s) that the haptic device should not exceed." );
//        readParameterInt(   &m_config.hapticsVelocityLimitPenaltyTime,                       "hapticsVelocityLimitPenaltyTime",                        1, "The duration of time that the haptic device does not get any non-zero force and torque commands after exceeding a velocity limit." );
//
//        readParameterBool(  &m_config.hapticsKukaMode,                                       "hapticsKukaMode",                                        1, "Whether we are currently using the plugin with a KUKA LBR robot." );
//        readParameterFloat(  m_config.hapticsKukaDrillHeadOffset,                            "hapticsKukaDrillHeadOffset",                             3, "The virtual tool mass that will be simualted by the virtual gravity (in kg)." );
//        readParameterFloat(  m_config.hapticsKukaCenterOfMass,                               "hapticsKukaCenterOfMass",                                3, "The virtual tool mass that will be simualted by the virtual gravity (in kg)." );
//        readParameterFloat( &m_config.hapticsKukaDrillMass,                                  "hapticsKukaDrillMass",                                   1, "The virtual tool mass that will be simualted by the virtual gravity (in kg)." );
//        readParameterBool(  &m_config.hapticsKukaUseVirtualGravity,                          "hapticsKukaUseVirtualGravity",                           1, "Whether to add a virtual gravity force to the rendered force." );
//        readParameterFloat(  m_config.hapticsKukaVirtualGravity,                             "hapticsKukaVirtualGravity",                              3, "The virtual tool mass that will be simualted by the virtual gravity (in kg)." );
//        readParameterFloat( &m_config.hapticsDampingLinear,                                  "hapticsDampingLinear",                                   1, "The linear damping of the KUKA robot in contact (as a fraction of the linear velocity). (0.2 means 20% of the robots velocity are converted to a force opposit of the velocity)" );
//        readParameterFloat( &m_config.hapticsDampingAngularOrigin,                           "hapticsDampingAngularOrigin",                            1, "Angular damping of the KUKA robot in contact (as a fraction of angular velocity at the robot flange). (High value => robot will try harder to keep the orientation of the flange the same)" );
//        readParameterFloat( &m_config.hapticsDampingAngularEndeffector,                      "hapticsDampingAngularEndeffector",                       1, "Angular damping of the KUKA robot in contact (as a fraction of angular velocity at the drill head). (High value => robot will try harder to keep the orientation of the drill head the same)\n# This easily leads to very high torque damping, as the drill head is at a large lever, consequently it naturally rotates a lot. (if enabled, choose very small value)" );
//        readParameterFloat( &m_config.hapticsDampingInFreeState,                             "hapticsDampingInFreeState",                              1, "The modifier that is used to reduce the effect of all damping when the robot is not in contact. (0.05 means 5% of the contact damping is applied when the robot is not in contact)\n# Could also be chosen to be 0." );
//
//        readParameterFloat( &m_config.transformationDeviceToChaiScale,                       "transformationDeviceToChaiScale",                        1, "How to scale the device positions of non-KUKA devices when transforming to Chai3D space." );
//        readParameterFloat( &m_config.transformationKukaToWorldScale,                        "transformationKukaToWorldScale",                         1, "How to scale the device positions of KUKA LBR robot when transforming to world space." );
//
//        readParameterBool(  &m_config.experimentMode,                                        "experimentMode",                                         1, "Whether to run the TUC experiment. (Only tested in Unreal)" );
//        readParameterFloat(  m_config.experimentToolToEnvironment,                           "experimentToolToEnvironment",                           16, "How fast the drill advances in the experiment (in m/s)." );
//        readParameterFloat(  m_config.experimentStartPushback,                               "experimentStartPushback",                                3, "How fast the drill advances in the experiment (in m/s)." );
//        readParameterFloat(  m_config.experimentDrillheadRadius,                             "experimentDrillheadRadius",                              3, "How fast the drill advances in the experiment (in m/s)." );
//        readParameterFloat( &m_config.experimentLeadVelocity,                                "experimentLeadVelocity",                                 1, "How fast the drill advances in the experiment (in m/s)." );
//
//        readParameterBool(  &m_config.nonlinearUseGrowingForce,                              "nonlinearUseGrowingForce",                               1, "Whether to use non-linear growing force (so basically, no matter the depth of penetration, the maximum force is never reached, however we guarantee that F(dop1) > F(dop2) if dop1 > dop2, we ideally it should never feel like the material gets softer after a certain point, since the force does not grow anymore (that's the case with linear growth)).\n# Does not work with hapticsForceMod." );
//        readParameterFloat( &m_config.nonlinearMaximumDepthOfPenetration,                    "nonlinearMaximumDepthOfPenetration",                     1, "The next three parameters define the non-linear force function.\n# The dop at which the maximum force should be reached." );
//        readParameterFloat( &m_config.nonlinearMaximumForceOfDevice,                         "nonlinearMaximumForceOfDevice",                          1, "The maximum force that the haptic device can render (adjust manually per device)." );
//        readParameterFloat( &m_config.nonlinearForceRatioAtMaxDepth,                         "nonlinearForceRatioAtMaxDepth",                          1, "The percent of the maximum force that should be rendered when reaching the earlier defined nonlinearMaximumDepthOfPenetration. (0.8 would mean 80% of the haptic devices maximum force will be rendered when the user penetrates as far as nonlinearMaximumDepthOfPenetration)\n# The function will asympotetly approach the maximum, but never actually reach it." );
//
//        readParameterBool(  &m_config.stictionUse,                                           "stictionUse",                                            1, "Whether to use stiction torque. (Stiction is when the drilling friction gets so high that the tool starts rotating because it gets stuck in the removable material, since it can't overcome the drilling friction)" );
//        readParameterFloat( &m_config.stictionDensityToStictionMod,                          "stictionDensityToStictionMod",                           1, "The percent of the density-weighted volume sum to convert to stiction. (High value => more stiction)" );
//        readParameterFloat( &m_config.stictionForceToStictionMod,                            "stictionForceToStictionMod",                             1, "The percent of the force magnitude to convert to stiction. (High value => more stiction)" );
//        readParameterFloat( &m_config.stictionMod,                                           "stictionMod",                                            1, "A general modifier to adjust the stiction strength. (High value => more stiction)" );
//
//        readParameterBool(  &m_config.hubbelUse,                                             "hubbelUse",                                              1, "Whether to use randomly occuring 'hubbels'" );
//        readParameterFloat( &m_config.hubbelDurationMean,                                    "hubbelDurationMean",                                     1, "Mean value of the randomly generated duration over which the 'hubbel' occures (in ms)" );
//        readParameterFloat( &m_config.hubbelDurationStdDev,                                  "hubbelDurationStdDev",                                   1, "Standard deviation of the randomly generated duration over which the 'hubbel' occures" );
//        readParameterFloat( &m_config.hubbelDensityMax,                                      "hubbelDensityMax",                                       1, "Standard deviation is derived from current contact surface. Low contact area => high std. deviation => high chance of large 'hubbel'\n# Minimal std. dev. determines how 'hubbel'-free good contact can be (lower means less 'hubbel')" );
//        readParameterFloat( &m_config.hubbelStdDevMin,                                       "hubbelStdDevMin",                                        1, "Maximum std. dev. determines how likely and hard the 'hubbels' will be in case of bad contact." );
//        readParameterFloat( &m_config.hubbelStdDevMax,                                       "hubbelStdDevMax",                                        1, "Determines what contact surface leads to minimal std. dev. for 'hubbel' generation. (Higher values increase 'hubbel' likelyhood & effect)" );
//        readParameterFloat( &m_config.hubbelIntensityMod,                                    "hubbelIntensityMod",                                     1, "Modifier for the intensity of hubbels" );
//        readParameterFloat( &m_config.hubbelMax,                                             "hubbelMax",                                              1, "Maximum 'hubbel' value that can be generated (this is used as multiplier for the force, hubbelMax 10 would mean a 'hubbel' could result in worst case in 10 times higher force than normal)" );
//        readParameterFloat( &m_config.hubbelProbabilityMod,                                  "hubbelProbabilityMod",                                   1, "Used to modify the likelyhood of a hubbel possibly occuring. (Higher value => higher probability)" );
//    }
//
//    Config* Config::getInstance()
//    {
//        static Config instance;
//        return &instance;
//    }
//
//    ConfigState &Config::Access()
//    {
//        return getInstance()->m_config;
//    }
//
//    void Config::readParameterBool(bool *targetVar, const char *name, int nrOfValues, const char *desc)
//    {
//        std::ifstream in( ConfigFileName );
//        std::string line;
//        bool found = false;
//
//        while (std::getline(in, line)) {
//            // Read current line
//            std::istringstream iss(line);
//            std::vector<std::string> tokens{ std::istream_iterator<std::string>{iss},
//                                             std::istream_iterator<std::string>{} };
//
//            // Skip empty lines
//            if (tokens.size() <= 0)
//                continue;
//
//            if (tokens[0] == name) {
//                found = true;
//                std::cout << "[Config] " << name << ":";
//                for (int i=0; i<nrOfValues; i++) {
//                    targetVar[i] = std::stoul(tokens[1+i]) != 0;
//                    std::cout << " " << (targetVar[i]?"1":"0");
//                }
//                std::cout << std::endl;
//                break;
//            }
//        }
//
//        in.close();
//
//        if (!found) {
//            std::cout << "[Config] (Not found => Writing default) " << name << ":";
//            std::ofstream out( ConfigFileName, std::ofstream::app );
//            out << "# " << desc << std::endl;
//            out << name;
//            for (int i=0; i<nrOfValues; i++) {
//                out << " " << targetVar[i];
//                std::cout << " " << targetVar[i];
//            }
//            out << std::endl;
//            std::cout << std::endl;
//            out.close();
//        }
//    }
//
//    void Config::readParameterInt(int *targetVar, const char *name, int nrOfValues, const char *desc)
//    {
//        std::ifstream in( ConfigFileName );
//        std::string line;
//        bool found = false;
//
//        while (std::getline(in, line)) {
//            // Read current line
//            std::istringstream iss(line);
//            std::vector<std::string> tokens{ std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{} };
//
//            // Skip empty lines
//            if (tokens.size() <= 0)
//                continue;
//
//            if (tokens[0] == name) {
//                found = true;
//                std::cout << "[Config] " << name << ":";
//                for (int i=0; i<nrOfValues; i++) {
//                    targetVar[i] = std::stoul(tokens[1+i]);
//                    std::cout << " " << targetVar[i];
//                }
//                std::cout << std::endl;
//                break;
//            }
//        }
//
//        in.close();
//
//        if (!found) {
//            std::cout << "[Config] (Not found => Writing default) " << name << ":";
//            std::ofstream out( ConfigFileName, std::ofstream::app );
//            out << "# " << desc << std::endl;
//            out << name;
//            for (int i=0; i<nrOfValues; i++) {
//                out << " " << targetVar[i];
//                std::cout << " " << targetVar[i];
//            }
//            out << std::endl;
//            std::cout << std::endl;
//            out.close();
//        }
//    }
//
//    void Config::readParameterFloat(float *targetVar, const char *name, int nrOfValues, const char *desc)
//    {
//        std::ifstream in( ConfigFileName );
//        std::string line;
//        bool found = false;
//
//        while (std::getline(in, line)) {
//            // Read current line
//            std::istringstream iss(line);
//            std::vector<std::string> tokens{ std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{} };
//
//            // Skip empty lines
//            if (tokens.size() <= 0)
//                continue;
//
//            if (tokens[0] == name) {
//                found = true;
//                std::cout << "[Config] " << name << ":";
//                for (int i=0; i<nrOfValues; i++) {
//                    targetVar[i] = std::stof(tokens[1+i]);
//                    std::cout << " " << targetVar[i];
//                }
//                std::cout << std::endl;
//                break;
//            }
//        }
//
//        in.close();
//
//        if (!found) {
//            std::cout << "[Config] (Not found => Writing default) " << name << ":";
//            std::ofstream out( ConfigFileName, std::ofstream::app );
//            out << "# " << desc << std::endl;
//            out << name;
//            for (int i=0; i<nrOfValues; i++) {
//                out << " " << targetVar[i];
//                std::cout << " " << targetVar[i];
//            }
//            out << std::endl;
//            std::cout << std::endl;
//            out.close();
//        }
//    }
//
//}
