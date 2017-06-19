#include <uWS/uWS.h>
#include <iostream>
#include "json.hpp"
#include <math.h>
#include "ukf.h"
#include "tools.h"


using namespace std;

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
std::string hasData(std::string s) 
{
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.find_first_of("]");
    if (found_null != std::string::npos) 
    {
        return "";
    }
    else if (b1 != std::string::npos && b2 != std::string::npos) 
    {
        return s.substr(b1, b2 - b1 + 1);
    }
    return "";
}

////////////////////////////////////////////////////////////////////////////////////////////////
// setup for simulator
////////////////////////////////////////////////////////////////////////////////////////////////
int main_forSimulator()
{
    uWS::Hub h;

    // Create a UKF instance
    UKF ukf;

    // used to compute the RMSE later
    Tools tools;
    vector<VectorXd> estimations;
    vector<VectorXd> ground_truth;

    // log values from simulator
    string out_file_name("../simulator_values.txt");
    ofstream output_file(out_file_name.c_str(), ofstream::out);
    if (!output_file.is_open())
    {
        cerr << "Cannot open output file: " << out_file_name << endl;
        exit(EXIT_FAILURE);
    }

    h.onMessage(
        [&ukf, &tools, &estimations, &ground_truth, &output_file]
        (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) 
    {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event

        if (length && length > 2 && data[0] == '4' && data[1] == '2')
        {

            auto s = hasData(std::string(data));
            if (s != "") 
            {
                auto j = json::parse(s);

                std::string event = j[0].get<std::string>();

                if (event == "telemetry") 
                {
                    // j[1] is the data JSON object

                    string sensor_measurment = j[1]["sensor_measurement"];

                    output_file << sensor_measurment << std::endl;

                    MeasurementPackage meas_package;
                    istringstream iss(sensor_measurment);
                    long long timestamp;

                    // reads first element from the current line
                    string sensor_type;
                    iss >> sensor_type;

                    if (sensor_type.compare("L") == 0) 
                    {
                        meas_package.sensor_type_ = MeasurementPackage::LASER;
                        meas_package.raw_measurements_ = VectorXd(2);
                        float px;
                        float py;
                        iss >> px;
                        iss >> py;
                        meas_package.raw_measurements_ << px, py;
                        iss >> timestamp;
                        meas_package.timestamp_ = timestamp;
                    }
                    else if (sensor_type.compare("R") == 0) 
                    {
                        meas_package.sensor_type_ = MeasurementPackage::RADAR;
                        meas_package.raw_measurements_ = VectorXd(3);
                        float ro;
                        float theta;
                        float ro_dot;
                        iss >> ro;
                        iss >> theta;
                        iss >> ro_dot;
                        meas_package.raw_measurements_ << ro, theta, ro_dot;
                        iss >> timestamp;
                        meas_package.timestamp_ = timestamp;
                    }
                    float x_gt;
                    float y_gt;
                    float vx_gt;
                    float vy_gt;
                    iss >> x_gt;
                    iss >> y_gt;
                    iss >> vx_gt;
                    iss >> vy_gt;
                    VectorXd gt_values(4);
                    gt_values(0) = x_gt;
                    gt_values(1) = y_gt;
                    gt_values(2) = vx_gt;
                    gt_values(3) = vy_gt;
                    ground_truth.push_back(gt_values);

                    //Call ProcessMeasurment(meas_package) for Kalman filter
                    ukf.ProcessMeasurement(meas_package);

                    //Push the current estimated x,y positon from the Kalman filter's state vector

                    VectorXd estimate(4);

                    VectorXd x = ukf.GetX();
                    double p_x = x(0);
                    double p_y = x(1);
                    double v = x(2);
                    double yaw = x(3);

                    double v1 = cos(yaw)*v;
                    double v2 = sin(yaw)*v;

                    estimate(0) = p_x;
                    estimate(1) = p_y;
                    estimate(2) = v1;
                    estimate(3) = v2;

                    estimations.push_back(estimate);

                    VectorXd RMSE = tools.CalculateRMSE(estimations, ground_truth);

                    json msgJson;
                    msgJson["estimate_x"] = p_x;
                    msgJson["estimate_y"] = p_y;
                    msgJson["rmse_x"] = RMSE(0);
                    msgJson["rmse_y"] = RMSE(1);
                    msgJson["rmse_vx"] = RMSE(2);
                    msgJson["rmse_vy"] = RMSE(3);
                    auto msg = "42[\"estimate_marker\"," + msgJson.dump() + "]";
                    // std::cout << msg << std::endl;
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                }
            }
            else 
            {
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
    });

    // We don't need this since we're not using HTTP but if it's removed the program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) 
    {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1)
        {
            res->end(s.data(), s.length());
        }
        else
        {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) 
    {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) 
    {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port))
    {
        std::cout << "Listening to port " << port << std::endl;
    }
    else
    {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();

    if (output_file.is_open())
    {
        output_file.close();
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// setup for data file
////////////////////////////////////////////////////////////////////////////////////////////////
void check_arguments(int argc, char* argv[])
{
    string usage_instructions = "Usage instructions: ";
    usage_instructions += argv[0];
    usage_instructions += " path/to/input.txt output.txt";

    bool has_valid_args = false;

    // make sure the user has provided input and output files
    if (argc == 1)
    {
        cerr << usage_instructions << endl;
    }
    else if (argc == 2)
    {
        cerr << "Please include an output file.\n" << usage_instructions << endl;
    }
    else if (argc == 3)
    {
        has_valid_args = true;
    }
    else if (argc > 3)
    {
        cerr << "Too many arguments.\n" << usage_instructions << endl;
    }

    if (!has_valid_args)
    {
        exit(EXIT_FAILURE);
    }
}

void check_files(
    ifstream& in_file, string& in_name,
    ofstream& out_file, string& out_name)
{
    if (!in_file.is_open())
    {
        cerr << "Cannot open input file: " << in_name << endl;
        exit(EXIT_FAILURE);
    }

    if (!out_file.is_open())
    {
        cerr << "Cannot open output file: " << out_name << endl;
        exit(EXIT_FAILURE);
    }
}

void readInputFile(
    ifstream& in_file_,
    vector<MeasurementPackage>& measurement_pack_list)
{
    string line;

    // prep the measurement packages (each line represents a measurement at a
    // timestamp)
    while (getline(in_file_, line))
    {
        string sensor_type;
        MeasurementPackage meas_package;
        istringstream iss(line);
        long long timestamp;

        // reads first element from the current line
        iss >> sensor_type;
        if (sensor_type.compare("L") == 0)
        {
            // LASER MEASUREMENT

            // read measurements at this timestamp
            meas_package.sensor_type_ = MeasurementPackage::LASER;
            meas_package.raw_measurements_ = VectorXd(2);
            float x;
            float y;
            iss >> x;
            iss >> y;
            meas_package.raw_measurements_ << x, y;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
            measurement_pack_list.push_back(meas_package);
        }
        else if (sensor_type.compare("R") == 0)
        {
            // RADAR MEASUREMENT

            // read measurements at this timestamp
            meas_package.sensor_type_ = MeasurementPackage::RADAR;
            meas_package.raw_measurements_ = VectorXd(3);
            float ro;
            float phi;
            float ro_dot;
            iss >> ro;
            iss >> phi;
            iss >> ro_dot;
            meas_package.raw_measurements_ << ro, phi, ro_dot;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
            measurement_pack_list.push_back(meas_package);
        }
    }
}

//int main(int argc, char* argv[])
int main_forDataFile(int argc, char* argv[])
{
    check_arguments(argc, argv);

    string in_file_name_ = argv[1];
    ifstream in_file_(in_file_name_.c_str(), ifstream::in);

    string out_file_name_ = argv[2];
    ofstream out_file_(out_file_name_.c_str(), ofstream::out);

    check_files(in_file_, in_file_name_, out_file_, out_file_name_);

    vector<MeasurementPackage> measurement_pack_list;

    readInputFile(in_file_, measurement_pack_list);

    // create UKF instance    
    UKF ukf;

    // used to compute the RMSE later
    vector<VectorXd> estimations;
    vector<VectorXd> ground_truth;

    //Call the EKF-based fusion
    size_t N = measurement_pack_list.size();
    for (size_t k = 0; k < N; ++k)
    {
        // start filtering from the second frame (the speed is unknown in the first frame)
        ukf.ProcessMeasurement(measurement_pack_list[k]);

        // output the estimation
        VectorXd x = ukf.GetX();
        out_file_ << x(0) << "\t";
        out_file_ << x(1) << "\t";
        out_file_ << x(2) << "\t";
        out_file_ << x(3) << "\t";

        // output the measurements
        if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER)
        {
            // output the estimation
            out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";
            out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";
        }
        else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR)
        {
            // output the estimation in the cartesian coordinates
            double ro = measurement_pack_list[k].raw_measurements_(0);
            double phi = measurement_pack_list[k].raw_measurements_(1);
            out_file_ << ro * cos(phi) << "\t"; // p1_meas
            out_file_ << ro * sin(phi) << "\t"; // ps_meas
        }

        estimations.push_back(x);
    }

    // close files
    if (out_file_.is_open())
    {
        out_file_.close();
    }

    if (in_file_.is_open())
    {
        in_file_.close();
    }

    ukf.CalculateNisConsistency();

    return 0;
}

int main(int argc, char* argv[])
{
    //main_forDataFile(argc, argv);
    main_forSimulator();
}
