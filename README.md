# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

## Introduction
The goal of this project is to drive a car in a highway scenario. In this project we interact with a simulator that includes our own vehicle, the highway and some traffic. The simulator provides us with telemetry data of all the vehicles and the highway waypoints.
During the driving we must comply with some requirements: 

1. The vehicle must drive at least 4.32 miles without incident.
2. The maximum speed must not be exceeded.
3. The maximum acceleration and jerk must not be exceeded.
4. The vehicle must not collide with other vehicles.
5. The vehicle must be able to change lane.
6. The vehicle must stay on its lane.

## Implementation

To achieve the goal of this project, using the provided waypoints and the telemetry data, I follow this steps:

### 1. Check for possible collisions

For this purpose, we read through all the telemetry data. For each vehicle, using the current position and speed, we predict its future position in the next two seconds. Then if at any intermediate steps the vehicles are closer than a threshold, a collision flag is raised. As a result, we get a vector that tell us if there is a potential collision on each lane. The method to check for collisions is checkCollision and can be found on line 289 of main.cpp.

### 2. Check for further traffic

To get in the fastest lane well ahead of encountering the vehicles, we check for traffic ahead. To do this, using the telemetry data, we check for the closest vehicle ahead on each lane. Then we predict the future distance in 5 seconds, this way even if the fastest vehicle is closer, we will be able to select the fastest path. This code is in lines 501 to 522 of main.cpp

### 3. Speed Control

In the approach I chose, for simplicity, I keep the speed and lane control apart. For the speed control, we check the closest vehicle on the current or target lane. If this vehicle is closer than 40 meters, we start decreasing the speed, if it's closer than 15 meter, we match the speed. Finally, if the vehicle is closer than 10 meter, we will further reduce the speed to increase the distance. This code is in lines 524 to 594 in main.cpp

### 4. Select Target Lane

Based on the current lane, collision and further traffic information, we calculate the cost for all the lanes, the lane with the lowest cost is set as the target lane. Basically we want to avoid all collisions, then we prioritize the lane with no traffic or the fastest traffic. While being on the side lanes, in case of a double lane change, we take into account collisions with the central lane.

Additionally, I implemented a counter to avoid too many lane changes, the vehicle is allowed to perform a lane change every four seconds. Also double lane changes are not allowed, these are made in two steps. Overall, this way of selecting the target lane is quite simple, easy to understand, and effective.

Below some examples of the car behaviour can be found.

In the first example we can see how the vehicle handles low traffic quite well, changing lane far before encountering the vehicles.

![Alt Text](https://github.com/cbg15589/CarND-Path-Planning-Project/blob/master/media/Low%20Traffic_fast.gif)

In the second example we see a heavier traffic situation, where the vehicle waits matching the front vehicle's speed until there is a safe window to overtake.

![Alt Text](https://github.com/cbg15589/CarND-Path-Planning-Project/blob/master/media/Low%20Traffic_2_fast.gif)

In the last example we can see a similar example, but this time after starting the overtake it doubts and goes back to lane 2. This is due to the car on lane 1 constantly braking and accelerating. This could be improved, but at least the overtake is done in a safe manner.

![Alt Text](https://github.com/cbg15589/CarND-Path-Planning-Project/blob/master/media/Heavy%20Traffic_fast.gif)

### 5. Trajectory Generation
   Finally, once we know on which lane we want to drive, it's time to generate the trajectory. Here I mainly used the project's walkthrough code with a minor tweak. Here we use the last point of the previous path and some of the highway waypoints further ahead to fit a spline. Based on the starting code, I had to dynamically adapt the distance between the waypoints based on the desired speed, without this, at low speeds, the vehicle didn't follow the center of the lane accurately.
   
The vehicle in the simulator travels to the next point every 0.02 so we need to split the spline into points with the required distance between them so that we achieved the demand speed. Then we add this points to the previous path and send the first 50 point to the simulator, using the previous path helps us to achieve smooth transitions.

## Conclusion

Overall, I'm quite satisfied with the results, achieving a personal best of 158.51 miles in 3:18:41 without any incident. This averages a speed of 47.87 mph, only two miles below the highway maximum speed.

![Alt Text](https://github.com/cbg15589/CarND-Path-Planning-Project/blob/master/media/Personal_Best.PNG)

Even then, all the incidents I found during a 12 hour run could not really be avoided without complex solutions, the incidents could be classified into 3 types:

#### 1.  Traffic Vehicle suddenly dissapearing from telemetry.

This could be avoided by introducing some data plausability checks, if a car suddenly disappears from the telemetry, we could at least predict a position based on previous observations.

![Alt Text](https://github.com/cbg15589/CarND-Path-Planning-Project/blob/master/media/Lost_Telemetry_Explained.PNG)

#### 2. Another vehicle changing to the same lane.

In this case, this incident could be avoided with emergency braking or aborting the lane change, either of those would result into an incident anyway.

![Alt Text](https://github.com/cbg15589/CarND-Path-Planning-Project/blob/master/media/Collision.gif)

#### 3. Lost connection with the simulator 

This happens very rarely, but sadly this is what happened on my personal best, suddenly the path planner lost connection with the simulator and the vehicle stopped.

Personal Best Video:
[Link](https://www.youtube.com/watch?v=5M1MSSQxhYk&feature=youtu.be)
   
   
