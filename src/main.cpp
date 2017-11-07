#define _USE_MATH_DEFINES

#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

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
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
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

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
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
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
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
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;
	while(s > maps_s[min(prev_wp+1, (int)(maps_s.size() - 1))] && (prev_wp < (int)(maps_s.size()-1)))
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

// Calculate trajectory
vector<vector<double>> generate_trajectory(double car_x, double car_y, double car_s, double car_yaw, double ref_vel, int lane, vector<double>previous_path_x, vector<double>previous_path_y, vector<double>map_waypoints_s, vector<double>map_waypoints_x, vector<double>map_waypoints_y) {
	vector<double> ptsx;
	vector<double> ptsy;

	// Reference x,y, yaw states
	// Either we will reference the starting point as where the car is or at the previous paths end point

	double ref_x = car_x;
	double ref_y = car_y;
	double ref_yaw = deg2rad(car_yaw);
	int prev_size = previous_path_x.size();

	if (prev_size<2)
	{
		// Use two points that make the path tangent to the car
		double prev_car_x = car_x - cos(car_yaw);
		double prev_car_y = car_y - sin(car_yaw);

		ptsx.push_back(prev_car_x);
		ptsx.push_back(car_x);

		ptsy.push_back(prev_car_y);
		ptsy.push_back(car_y);
	}
	// Use the previous path's end point as starting reference
	else
	{
		// redefine reference state as previous path end point
		ref_x = previous_path_x[prev_size - 1];
		ref_y = previous_path_y[prev_size - 1];

		double ref_x_prev = previous_path_x[prev_size - 2];
		double ref_y_prev = previous_path_y[prev_size - 2];
		ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

		// Use two points that make the path tangent to the previous path's end point
		ptsx.push_back(ref_x_prev);
		ptsx.push_back(ref_x);

		ptsy.push_back(ref_y_prev);
		ptsy.push_back(ref_y);

	}

	// In Frenet add evenly spaced points ahead of the starting reference,
	// space between points depends on the vehicle speed. Lower speeds requires
	// less distance between waypoints for an accurate trayectory that stays in lane.
	vector<double> next_wp1 = getXY(car_s + max((40*ref_vel/49.5),(double)10), (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
	vector<double> next_wp2 = getXY(car_s + max((80 * ref_vel / 49.5), (double)20), (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
	vector<double> next_wp3 = getXY(car_s + max((120 * ref_vel / 49.5), (double)30), (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

	ptsx.push_back(next_wp1[0]);
	ptsx.push_back(next_wp2[0]);
	ptsx.push_back(next_wp3[0]);

	ptsy.push_back(next_wp1[1]);
	ptsy.push_back(next_wp2[1]);
	ptsy.push_back(next_wp3[1]);

	for (int i = 0; i < ptsx.size(); i++)
	{
		// Shift car reference angle to 0 degrees
		double shift_x = ptsx[i] - ref_x;
		double shift_y = ptsy[i] - ref_y;

		ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
		ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
	}
	// Create a spline
	tk::spline s;

	// Set (x,y) points to the spline
	s.set_points(ptsx, ptsy);

	// Define the actual (x,y) points we will use for planner
	vector<double> next_x_vals;
	vector<double> next_y_vals;

	// Start with all the previous path points from last time
	for (int i = 0; i < previous_path_x.size(); i++)
	{
		next_x_vals.push_back(previous_path_x[i]);
		next_y_vals.push_back(previous_path_y[i]);
	}

	// Calculate how to break up spline so that we travel at our desired reference velocity
	double target_x = 30.0;
	double target_y = s(target_x);
	double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));

	double x_add_on = 0;

	// Fill up the rest of our path planner after filling it with previous points, here we will always output 50 points

	for (int i = 1; i <= 50 - next_x_vals.size(); i++)
	{
		double N = (target_dist / (0.02*ref_vel / 2.24));
		double x_point = x_add_on + (target_x) / N;
		double y_point = s(x_point);

		x_add_on = x_point;

		double x_ref = x_point;
		double y_ref = y_point;

		// Rotate back to normal after rotating it earlier
		x_point = (x_ref * cos(ref_yaw) - y_ref*sin(ref_yaw));
		y_point = (x_ref * sin(ref_yaw) + y_ref*cos(ref_yaw));

		x_point += ref_x;
		y_point += ref_y;

		next_x_vals.push_back(x_point);
		next_y_vals.push_back(y_point);
	}
	vector<vector<double>> trajectory;
	trajectory.push_back(next_x_vals);
	trajectory.push_back(next_y_vals);

	return trajectory;

}

// Check Collision with other vehicle
bool checkCollision(double s1, double s2, double v1, double v2) {
	
	bool collision = false;
	// Time horizon to look for collision
	double time_gap = 2;
	// Predict vehicles positions between present time and the time horizon
	for (int a = 0; a < 101; a++)
	{
		double s1_next = s1 + v1*time_gap*(a/100);
		double s2_next = s2 + v2*time_gap*(a/100);
		// If the distance between vehicles is smaller that 5 at time = 0,
		// or 15 at time = 2. A possible collision is flaged
		if (abs(s1_next - s2_next) < (5 + 10*(a/100)))
		{
			collision = true;
		}
	}

	return collision;
}

// Check if there id traffic far ahead (Not used)
bool checkFarTraffic(double s1, double s2, double v1, double v2, double time_gap) {

	bool traffic = false;
	double s1_next = s1 + v1*time_gap;
	double s2_next = s2 + v2*time_gap;
	if (s1_next > s2_next && s1 < s2)
	{
		traffic = true;
	}

	return traffic;
}

// Cost calculation used to select the target lane.
vector<double> cost_calc(int lane, vector<bool>laneCollision, vector<double>closest_car_distance) {
		vector<double> cost = { 0,0,0 };
		double collision_cost = 10000;

		// Cost assigned for traffic far ahead
		// Distances are sorted, the lowest cost is assignt to the lane with the furthest car
		// In case of no traffic, this will favour right lane driving.
		vector<double> lane_traffic_cost = { 0,0,0 };
		vector<std::pair<double, int>> sorted_distances = { { closest_car_distance[0], 0 },{ closest_car_distance[1], 1 },{ closest_car_distance[2], 2 } };
		sort(sorted_distances.begin(), sorted_distances.end());
		lane_traffic_cost[sorted_distances[0].second] = 1000;
		lane_traffic_cost[sorted_distances[1].second] = 500;
		lane_traffic_cost[sorted_distances[2].second] = 0;

		double lane_change_right_cost = 5;
		double lane_change_left_cost = 10;
	
		if (lane == 0) {
			cost[0] = min(collision_cost, laneCollision[0] * collision_cost + lane_traffic_cost[0]);
			cost[1] = lane_change_right_cost + min(collision_cost, laneCollision[1] * collision_cost + lane_traffic_cost[1]);
			cost[2] = lane_change_right_cost * 2  + laneCollision[1] * collision_cost + laneCollision[2] * collision_cost + lane_traffic_cost[2];
		}
		else if (lane == 1) {
			cost[0] = lane_change_left_cost + min(collision_cost, laneCollision[0] * collision_cost + lane_traffic_cost[0]);
			cost[1] = min(collision_cost, laneCollision[1] * collision_cost + lane_traffic_cost[1]);
			cost[2] = lane_change_right_cost + min(collision_cost, laneCollision[2] * collision_cost + lane_traffic_cost[2]);
		}
		else{
			cost[0] = lane_change_left_cost * 2 + laneCollision[1] * collision_cost + laneCollision[0] * collision_cost + lane_traffic_cost[0];;
			cost[1] = lane_change_left_cost + min(collision_cost, laneCollision[1] * collision_cost + lane_traffic_cost[1]);
			cost[2] = min(collision_cost, laneCollision[2] * collision_cost + lane_traffic_cost[2]);
		}
		return cost;
}

// Get lane based on the distance from the center of the road
double get_current_lane(double car_d) {
	int currentLane;
	if (car_d <= 4) currentLane = 0;
	if (car_d > 4 && car_d <= 8) currentLane = 1;
	if (car_d > 8) currentLane = 2;

	return currentLane;
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

  //start in lane 1
  int lane = 1;
  // Have a reference velocity to target
  double ref_vel = 0;
  // Counter for lane change
  int counter = 0;
  bool lane_changed = false;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &lane, &ref_vel, &counter,&lane_changed](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
			car_speed = car_speed / 2.24;// to m/s
			// Car s coordinate used for collision and trafic calculations
			double original_car_s = car_s;

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];


			// Set car s coordinate as the last value of the previous path. 
			// This is used for trayectory calculation.
			int prev_size = previous_path_x.size();

			if (prev_size > 0) {
				car_s = end_path_s;
			}

			// Flags used for speed control
			bool too_close = false;			
			bool match_speed = false;
			bool reduce_speed = false;

			// Variables used for lane change control
			vector<bool> laneCollision = { false, false, false };
			vector<double> closest_car_distance = { 9999, 9999, 9999 };
			vector<double> closest_car_distance_5s = { 9999, 9999, 9999 };
			vector<double> closest_car_speed = { 0, 0, 0 };

			// Speed and lane of the car in front.
			double front_car_speed = 0;
			int currentLane = get_current_lane(car_d);

			// Loop through all the telemetry data
			for (int i = 0; i < sensor_fusion.size(); i++)
			{
				double vx = sensor_fusion[i][3];
				double vy = sensor_fusion[i][4];
				double check_speed = sqrt(vx*vx + vy*vy);
				double check_car_s = sensor_fusion[i][5];
				float d = sensor_fusion[i][6];
				bool collision = false;

				// Check if there is a probable collision with this car and raise flag
				collision = checkCollision(original_car_s, check_car_s, car_speed, check_speed);

				if (collision == true)
				{
					if (d <= 4) laneCollision[0] = true;
					if (d > 4 && d <= 8) laneCollision[1] = true;
					if (d > 8) laneCollision[2] = true;
				}

				// If the car is ahead check distance, store closest car distance and speed in each lane.
				// Distance is calculated 5 seconds in the future to take into account the speed.
				double distance = (check_car_s - original_car_s);
				if (distance > 0 && distance < 150)
				{
					if (d <= 4) {
						if (closest_car_distance[0] > distance) {
							closest_car_distance[0] = distance;
							closest_car_distance_5s[0] = distance + check_speed * 5;
							closest_car_speed[0] = check_speed;
						}
					}
					else if (d > 4 && d <= 8) {
						if (closest_car_distance[1] > distance) {
							closest_car_distance[1] = distance;
							closest_car_distance_5s[1] = distance + check_speed * 5;
							closest_car_speed[1] = check_speed;
						}
					}
					else {
						if (closest_car_distance[2] > distance) {
							closest_car_distance[2] = distance;
							closest_car_distance_5s[2] = distance + check_speed * 5;
							closest_car_speed[2] = check_speed;
						}
					}
				}

				// If the car is ahead, close and on the current or target lane, regulate spees
				if ((d < (2 + 4 * currentLane + 2) && d >(2 + 4 * currentLane - 2)) || (d < (2 + 4 * lane + 2) && d >(2 + 4 * lane - 2)))
				{
					// Checks if car is ahead
					if (check_car_s > original_car_s)
					{
						// If car is closer than 40m, start reducing the speed
						if ((check_car_s - original_car_s) < 40)
						{
							front_car_speed = check_speed;
							reduce_speed = true;
						}
						// If car is 15m ahead match it's speed
						if ((check_car_s - original_car_s) < 15)
						{
							match_speed = true;
							laneCollision[currentLane] = true;
						}
						// If car is 10m ahead further reduce the speed
						if ((check_car_s - original_car_s) < 10)
						{
							too_close = true;
						}
						
					}
					
				}

			}
			// Report speeds for debugging
			cout << "Car Speed: " << car_speed << endl;
			cout << "Closest Car Speed: " << front_car_speed << endl;
			
			// Speed Control based on flags
			if (reduce_speed) {
				if (match_speed){
					if (too_close)
					{
						if ((ref_vel > front_car_speed*2.24*0.9))
						{
							ref_vel -= .224;
							cout << "Braking, car is closer than 10m" << endl;
						}
						else if (ref_vel < 49.5)
						{
							ref_vel += .224;
						}
					}
					else if (ref_vel > front_car_speed*2.24)
					{
						cout << "Matching speed with front car" << endl;
						ref_vel -= .224;
					}
					else if (ref_vel < 49.5) {
						ref_vel += .224;
					}
				}
				else if (ref_vel > front_car_speed*2.24*1.05)
				{
					cout << "Gradually reduce speed, car is closer that 40m" << endl;
					ref_vel -= .224;
				}
				else if (ref_vel < 49.5)
				{
					ref_vel += .224;
				}
				
			} else if (ref_vel < 49.5)
			{
				ref_vel += .224;
			}

			// If no previous lane change has been done, calculate costs for potential changes 
			if (lane_changed == false)
			{
				vector<double> cost = cost_calc(lane, laneCollision, closest_car_distance_5s);

				// Report status for debugging
				cout << "Current Lane: " << currentLane << "    Target Lane: " << lane << " Collissions:  " << laneCollision[0] << "  " << laneCollision[1] << "  " << laneCollision[2] << " Distance:  " << closest_car_distance[0] << "  " << closest_car_distance[1] << "  " << closest_car_distance[2] << " Distance_5sec:  " << closest_car_distance_5s[0]<< "  " << closest_car_distance_5s[1] << "  " << closest_car_distance_5s[2] << " Cost:  " << cost[0] << " " << cost[1] << " " << cost[2] << endl;
				
				// Set target lane as the one with the lowest cost
				int new_lane = min_element(cost.begin(), cost.end()) - cost.begin();

				// If lane change is required, flag it. Double lane changes are done in two steps
				// for confort.
				if (new_lane != lane) {
					lane_changed = true;
					cout << "Lane Change!" << endl;
					if (abs(lane - new_lane) == 2) {
						lane = 1;
					}
					else {
						lane = new_lane;
					}
					
				}
			} else {
				// Once a lane change is madde a counter is activated to prevent further lane
				// changes in the next 4 seconds (200*0.02)
				counter += 1;
				cout << counter << endl;
				if (counter >= 200) {
					lane_changed = false;
					counter = 0;
				}
			}

			// Generate new trajectory based on new target lane.
			vector<vector<double>> trajectory = generate_trajectory(car_x, car_y, car_s, car_yaw, ref_vel, lane, previous_path_x, previous_path_y, map_waypoints_s, map_waypoints_x, map_waypoints_y);
			
			// Send datato the simulator
          	json msgJson;

          	msgJson["next_x"] = trajectory[0];
          	msgJson["next_y"] = trajectory[1];

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

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
  if (h.listen("127.0.0.1",port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
