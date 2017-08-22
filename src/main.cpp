#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

#define LARGE_NUM 100000
#define MAX_SPEED 22.0
#define MAX_PATH_SIZE 50.0

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = LARGE_NUM;
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;
  
  // FSM and planner state variables
  bool changing_lanes = false;
  double target_lane = -1;
  double target_speed = -1;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&target_speed, &target_lane, &changing_lanes, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
			
			// traffic state variables
			bool traffic_ahead = false;
			bool traffic_left_frt = false;
			bool traffic_left_bck = false;
			bool traffic_right_frt = false;
			bool traffic_right_bck = false;
			
			// prepare to change lanes bool flag
			bool prep_change_lanes = false;
			
			// distance and speed states for other cars in traffic
			double min_diff_ahead = LARGE_NUM;
			double min_diff_left_frt = LARGE_NUM;
			double min_diff_left_bck = -LARGE_NUM;
			double min_diff_right_frt = LARGE_NUM;
			double min_diff_right_bck = -LARGE_NUM;
			
			double speed_ahead = MAX_SPEED;
			double speed_left_frt = MAX_SPEED;
			double speed_left_bck = MAX_SPEED;
			double speed_right_frt = MAX_SPEED;
			double speed_right_bck = MAX_SPEED;
			
			// car state variables			
			auto current_lane = floor(car_d/4.0);
			auto current_speed = car_speed;
			auto prev_path_size = previous_path_x.size();
			auto ref_y = car_y;
			auto ref_x = car_x;
			auto ref_yaw = car_yaw;
			auto ref_s = car_s;
			auto prev_ref_x = car_x - cos(car_yaw);
			auto prev_ref_y = car_y - sin(car_yaw);
			
			// if in changing lanes state, but already reached target_lane, switch to not changing lanes state
			if(changing_lanes && current_lane == target_lane){
				changing_lanes = false;
				cout << endl;
			}
			
			// use previously calculated path for a smoother route
			if(prev_path_size >= 2){
				prev_ref_x = previous_path_x[prev_path_size - 2];
				prev_ref_y = previous_path_y[prev_path_size - 2];
				ref_x = previous_path_x[prev_path_size - 1];
				ref_y = previous_path_y[prev_path_size - 1];
				ref_yaw = atan2(ref_y - prev_ref_y, ref_x - prev_ref_x);
				ref_s = end_path_s;
				current_lane = floor(end_path_d/4.0);
				current_speed = distance((ref_x - prev_ref_x), (ref_y - prev_ref_y), 0, 0);
				current_speed = current_speed/0.02;
			}
			
			// min and max distances from other cars in traffic to guarantee a safer route
			double min_diff_bck = -7.0;
			double max_diff_bck = 15.0;
			double max_diff_frt = 30.0;
			double max_diff_frt_sides = 60.0;
			
			// go through sensor_fusion list and record distances for closest cars
			// set traffic state variables
			for(auto i = 0; i < sensor_fusion.size(); i++){
				vector<double> vehicle = sensor_fusion[i];
				auto vehicle_speed = distance(vehicle[3], vehicle[4], 0, 0);
				auto vehicle_lane = floor(vehicle[6]/4.0);
				auto vehicle_s = vehicle[5];
				
				auto vehicle_s_fwd = vehicle_s + prev_path_size*0.02*vehicle_speed;
				
				auto s_diff = vehicle_s_fwd - ref_s;
				
				// if detected vehicle in lane to the left
				if(vehicle_lane == current_lane - 1){
					// check if vehicle meets safety conditions
					if(s_diff > min_diff_bck && s_diff < max_diff_frt_sides){
						// check if vehicle is in front or besides/behind the car
						if(s_diff > max_diff_bck){
							// check if it's the closest vehicle in front to the left
							if(s_diff < min_diff_left_frt){
								min_diff_left_frt = s_diff;
								speed_left_frt = vehicle_speed;
								traffic_left_frt = true;
							}
						} else {
							// check if it's the closest vehicle besides/behind to the left
							if(s_diff > min_diff_left_bck){
								min_diff_left_bck = s_diff;
								speed_left_bck = vehicle_speed;
								traffic_left_bck = true;
							}
						}
					}
				}
				// if detected vehicle in lane to the right
				if(vehicle_lane == current_lane + 1){
					// check if vehicle meets safety conditions
					if(s_diff > min_diff_bck && s_diff < max_diff_frt_sides){
						// check if vehicle is in front or besides/behind the car
						if(s_diff > max_diff_bck){
							// check if it's the closest vehicle in front to the right
							if(s_diff < min_diff_right_frt){
								min_diff_right_frt = s_diff;
								speed_right_frt = vehicle_speed;
								traffic_right_frt = true;
							}
						} else {
							// check if it's the closest vehicle besides/behind to the right
							if(s_diff > min_diff_right_bck){
								min_diff_right_bck = s_diff;
								speed_right_bck = vehicle_speed;
								traffic_right_bck = true;
							}
						}
					}
				}
				// if detected vehicle in same lane
				if(vehicle_lane == current_lane){
					// check if vehicle meets safety conditions
					if(s_diff > 0 && s_diff < max_diff_frt){
						// check if it's the closest vehicle in front
						if(s_diff < min_diff_ahead){
							min_diff_ahead = s_diff;
							speed_ahead = vehicle_speed;
							traffic_ahead = true;
						}
					}
				}
			}
			
			// if not changing lanes, keep lane unchanged and 
			// change speed to reflect the most up to date traffic state
			if(!changing_lanes){
				target_lane = current_lane;
				target_speed = speed_ahead;
			}
			
			// when there's traffic ahead
			if(traffic_ahead){
				
				// safety conditions for closest distance allowed 
				// and desired distance to be from the vehicle in front
				double critical_diff_ahead = 7.0;
				double desired_diff_ahead = 20.0;
				// if car doesn't meet safety conditions, set target speed accordingly
				if(min_diff_ahead < critical_diff_ahead || min_diff_ahead > desired_diff_ahead){
					target_speed = min(MAX_SPEED, target_speed + (min_diff_ahead - desired_diff_ahead)/((MAX_PATH_SIZE - prev_path_size)*0.02));
				}
				
				// planner's brain - it checks several conditions for the traffic state variables and
				// sets the FSM state and the targets (lane and speed) accordingly
				if(current_lane > 0 && !traffic_left_bck && !changing_lanes){
					if(!traffic_left_frt){
						target_lane = current_lane - 1;
						target_speed = speed_left_frt;
						//cout << "left0"; //debug
						prep_change_lanes = true;
					} else {
						if(min_diff_left_frt > min_diff_ahead + critical_diff_ahead){
							target_lane = current_lane - 1;
							target_speed = speed_left_frt;
							target_speed = min(MAX_SPEED, target_speed + (min_diff_left_frt - critical_diff_ahead)/((MAX_PATH_SIZE - prev_path_size)*0.02));
							//cout << "left1"; //debug
							prep_change_lanes = true;
						} else {
							if(speed_left_frt > target_speed + 3){
								target_lane = current_lane - 1;
								target_speed = speed_left_frt;
								//cout << "left2"; //debug
								prep_change_lanes = true;
							}
						}
					}
				}
				if(current_lane < 2 && !traffic_right_bck && !changing_lanes){
					if(!traffic_right_frt){
						target_lane = current_lane + 1;
						target_speed = speed_right_frt;
						//cout << "right0"; //debug
						prep_change_lanes = true;
					} else {
						if(min_diff_right_frt > min_diff_ahead + critical_diff_ahead){
							if(prep_change_lanes){
								if(min_diff_right_frt > min_diff_left_frt + critical_diff_ahead){
									target_lane = current_lane + 1;
									target_speed = speed_right_frt;
									target_speed = min(MAX_SPEED, target_speed + (min_diff_right_frt - critical_diff_ahead)/((MAX_PATH_SIZE - prev_path_size)*0.02));
									//cout << "right1"; //debug
								}
							} else {
								target_lane = current_lane + 1;
								target_speed = speed_right_frt;
								target_speed = min(MAX_SPEED, target_speed + (min_diff_right_frt - critical_diff_ahead)/((MAX_PATH_SIZE - prev_path_size)*0.02));
								//cout << "right2"; //debug
								prep_change_lanes = true;
							}
						} else {
							if(!prep_change_lanes){
								if(speed_right_frt > target_speed + 3){
									target_lane = current_lane + 1;
									target_speed = speed_right_frt;
									//cout << "right3"; //debug
									prep_change_lanes = true;
								}
							} else {
								if(min_diff_left_frt <= min_diff_ahead + critical_diff_ahead){
									if(speed_right_frt > speed_left_frt){
										target_lane = current_lane + 1;
										target_speed = speed_right_frt;
										//cout << "right4"; //debug
									}
								}
							}
						}
					}
				}
				if(current_lane == 1 && traffic_left_frt && traffic_right_frt && !traffic_left_bck && !traffic_right_bck && !changing_lanes && !prep_change_lanes){
					vector<double> traffic_speeds = {speed_ahead, speed_left_frt, speed_right_frt};
					sort(traffic_speeds.begin(), traffic_speeds.end());
					if(fabs(traffic_speeds[0] - target_speed) > 10){
						if(traffic_speeds[0] == speed_left_frt){
							target_lane = current_lane - 1;
							target_speed = speed_left_frt;
							//cout << "left3"; //debug
							prep_change_lanes = true;
						}
						if(traffic_speeds[0] == speed_right_frt){
							target_lane = current_lane + 1;
							target_speed = speed_right_frt;
							//cout << "right5"; //debug
							prep_change_lanes = true;
						}
					}
				}
			} else {
				if(current_lane == 0 && !traffic_right_bck && !changing_lanes){
					if(min_diff_right_frt > max_diff_frt){
						//cout << "right6"; //debug
						target_lane = current_lane + 1;
						prep_change_lanes = true; //make sure this is necessary!!! maybe safer without?
					}
				}
				if(current_lane == 2 && !traffic_left_bck && !changing_lanes){
					if(min_diff_left_frt > max_diff_frt){
						//cout << "left4"; //debug
						target_lane = current_lane - 1;
						prep_change_lanes = true; //make sure this is necessary!!! maybe safer without?
					}
				}
			}
			
			// if not changing lanes, check if planner is preparing to change
			if(!changing_lanes){
				changing_lanes = prep_change_lanes;
			}
			
			// use spline library to calculate smooth path given the targets lane and speed
			vector<double> spline_ptsx;
			vector<double> spline_ptsy;
			
			// use at the beginning of the spline the car state variables calculated before
			spline_ptsx.push_back(prev_ref_x);
			spline_ptsx.push_back(ref_x);
			spline_ptsy.push_back(prev_ref_y);
			spline_ptsy.push_back(ref_y);
			
			// set 2 more spline points at a horizon of 30m and 60m
			vector<double> next_spline_pt0 = getXY(ref_s + 30, (2 + 4*target_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_spline_pt1 = getXY(ref_s + 60, (2 + 4*target_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			
			spline_ptsx.push_back(next_spline_pt0[0]);
			spline_ptsx.push_back(next_spline_pt1[0]);
			
			spline_ptsy.push_back(next_spline_pt0[1]);
			spline_ptsy.push_back(next_spline_pt1[1]);
			
			// transform the points from global coordinates to car coordinates
			for(auto i = 0; i < spline_ptsx.size(); i++){
				double shift_x = spline_ptsx[i] - ref_x;
				double shift_y = spline_ptsy[i] - ref_y;
				
				spline_ptsx[i] = (shift_x*cos(-ref_yaw) - shift_y*sin(-ref_yaw));
				spline_ptsy[i] = (shift_x*sin(-ref_yaw) + shift_y*cos(-ref_yaw));
			}
			
			// calculate spline
			tk::spline s;
			s.set_points(spline_ptsx, spline_ptsy);

			// path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
			
			// use remaining coordinates from previously calculated path
			for(auto i = 0; i < prev_path_size; i++){
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			}
			
			// calculate step distances, given the target speed, the time step and
			// the safety limit for the acceleration
			double current_x = 0;
			for(auto i = 0; i < MAX_PATH_SIZE - prev_path_size; i++){
				if(current_speed < target_speed){
					current_speed = min(target_speed, current_speed + 0.15);
				} else {
					current_speed = max(target_speed, current_speed - 0.15);
				}
				double target_x = current_x + current_speed*0.02;
				double target_y = s(target_x);
				
				current_x = target_x;
				
				// transform points back from car coordinates to global coordinates
				double x_val = (target_x*cos(ref_yaw) - target_y*sin(ref_yaw));
				double y_val = (target_x*sin(ref_yaw) + target_y*cos(ref_yaw));
				
				x_val += ref_x;
				y_val += ref_y;
				
				next_x_vals.push_back(x_val);
				next_y_vals.push_back(y_val);
			}

			json msgJson;
          	// send path as json message to simulator
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































