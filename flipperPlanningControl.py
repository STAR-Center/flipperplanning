# -*- coding: utf-8 -*-
#!/usr/bin/env python2
from __future__ import division
import rospy
import numpy as np
from message_filters import Subscriber,ApproximateTimeSynchronizer
from sensor_msgs.msg import Imu
from geometry_msgs.msg import PoseStamped, Pose, Quaternion
from dynamixel_workbench_controllers.msg import XM2
from nav_msgs.msg import Path
from std_msgs.msg import Float64MultiArray, Header
from joy_control.msg import Float64StampedMultiArray
import tf
import copy
#from traverse_d import CSpace
from cspace import CSpace 

#from apscheduler.scheduler import Scheduler

#sched = Scheduler()
#sched.start()


first_frame = True

first_xm3 = 0
first_xm4 = 0
cur_xm3 = 0
cur_xm4 = 0
cur_wheel_pose = PoseStamped()

first_posi = None
first_euler = None
cur_vrpn_pose = PoseStamped()
CAR_WIDTH = 0.12

# the xm3, xm4 when we plan and run
st_xm3 = 0
st_xm4 = 0
#in_planning = False # if we are in planning run
on_planning_run = False# first frame we on planning run

#path_lst = []# path that will be updated in the first trigger frame
current_pt = None
next_pt = None
last_point = -1 #-1 next is not last, 0 next is last, 1 curr is last
last_m = 0 #m from start to currentpoint
first_gt_R = True # first time a > R

D_START = 0.40#0.40
H_START = 0.095-0.049#0.028#0.095 - 0.0
R = 0.035

REACH_ESPILON = 1e-2
STEADY_VELOCITY = -0.1
STEADY_ANGULAR = 0#angular = 0 means not rotate
DELTA_ANGLE_CODE = 0.01

alpha_commandCode = 0.
back_alpha_commandCode = 0.
velocity = 0.
angular = 0.

global control_pub, record_param_pub


def alpha2commandCode(alpha_in):
    alpha = min(56/180*np.pi, alpha_in)
    alpha = max(-np.pi/2, alpha)
    code = (-alpha+56/180*np.pi) / np.pi
    code = code if code < 0.81111 else 0.81111
    return code
def commandCode2alpha(comm):
    alpha = 56/180*np.pi - comm * np.pi
    return alpha
def reach_target(delta_m, wanted_delta_m):
	if delta_m - wanted_delta_m > 1e-6:
		return True
	else:
		return False

def publish():
    global velocity, angular, alpha_commandCode, back_alpha_commandCode, control_pub
    command = Float64MultiArray()
    rospy.loginfo("command is %f, %f"%(alpha_commandCode, back_alpha_commandCode))
    command.data.append(velocity)
    command.data.append(angular) 
    '''
    command.data.append(alpha_commandCode)#(front_l)#note that the range is 0~1 that is assume to be pi/2~-pi/2
    command.data.append(alpha_commandCode)#(front_r) 
    command.data.append(back_alpha_commandCode)#(back_l) 
    command.data.append(back_alpha_commandCode)#(back_r) 
    '''
    command.data.append(back_alpha_commandCode)#(front_l)#note that the range is 0~1 that is assume to be pi/2~-pi/2
    command.data.append(back_alpha_commandCode)#(front_r) 
    command.data.append(alpha_commandCode)#(back_l) 
    command.data.append(alpha_commandCode)#(back_r)  

    control_pub.publish(command) 

def publish_param():
    '''
        follow Soeren's request, we should also publish the target tuple, current alpha, beta.
        array[6+2]
    '''
    global next_pt, alpha_commandCode, back_alpha_commandCode
    record_data = Float64StampedMultiArray()
    record_data.header = Header()
    record_data.header.stamp = rospy.get_rostime()
    record_data.header.frame_id = 'robotInnerConfig'
    for i in range(6):
        record_data.array.data.append(next_pt[i])
    record_data.array.data.append(commandCode2alpha(alpha_commandCode))
    record_data.array.data.append(commandCode2alpha(back_alpha_commandCode))
    record_param_pub.publish(record_data)

def callback(xm3, xm4):
    global st_xm3, st_xm4, in_planning, on_planning_run, path_lst, control_pub, current_pt, next_pt, alpha_commandCode, back_alpha_commandCode, velocity, angular, last_point, first_gt_R, last_m#, sched
    #velocity = STEADY_VELOCITY
    #angular = STEADY_ANGULAR 


    if in_planning:
        if not on_planning_run:
            # in the first frame we make plan and get path
            on_planning_run = True
            st_xm3 = xm3.Position_Trajectory
            st_xm4 = xm4.Position_Trajectory

            current_pt = (D_START, R, 56/180*np.pi, 56*180/np.pi, 0., 0.145)#path_lst.pop()
            next_pt = path_lst.pop()

            last_m = 0

            d,a,alpha,back_alpha,theta,t = current_pt

            # since the upper bound is 56, so we the allowed code is in 56 to -90
            # the inital angle code should be 0
            alpha_commandCode = alpha2commandCode(alpha)
            back_alpha_commandCode = alpha2commandCode(back_alpha)
            #sched.add_interval_job(publish, seconds = 0.01)
        else:
            m_right = float(xm3.Position_Trajectory - st_xm3) / (-210) * 0.01
            m_left = float(xm4.Position_Trajectory - st_xm4) / (208) * 0.01

            m_right = (m_right+m_left)/2

            #rospy.loginfo("m %f"%m_left)
            d,a,alpha,back_alpha,theta, t = next_pt
            # compute the m_minus
            if a == R:
                m_minus = -(next_pt[0]-current_pt[0]) + (next_pt[4] - current_pt[4])*R
                rospy.loginfo("m_mimus a==r %f, d0 %f, d1 %f"%(m_minus,current_pt[0], next_pt[0]))
            else:
                m_minus = -(next_pt[5] - current_pt[5])
                rospy.loginfo("m_mimus a>r %f"%m_minus)

            '''
            if m_right-(D_START-current_pt[0])>0:
                velocity=0.
                angular = 0.
                current_pt = next_pt
                if len(path_lst) == 0:
                    if last_point == -1:
                        last_point=0
                    elif last_point == 0:
                        last_point = 1
                else:
                    next_pt = path_lst.pop()
            '''
            
            if reach_target(m_right-last_m, m_minus):
                rospy.loginfo('########################################################################################m_right %f, laft_m %f, m_minus %f'%(m_right, last_m, m_minus))
                last_m += m_minus
                current_pt = next_pt
                if len(path_lst) == 0:
                    if last_point == -1:
                        last_point=0
                    elif last_point == 0:
                        last_point = 1
                else:
                    next_pt = path_lst.pop()
                if last_point == 1:
                    velocity = 0.
                #angular = 0.
                    #alpha_commandCode=0
                    #back_alpha_commandCode=0
            else:
                velocity = STEADY_VELOCITY
                angular = STEADY_ANGULAR
            d,a,alpha,back_alpha,theta, t = next_pt
            next_alpha_commandCode = alpha2commandCode(alpha) 
            rospy.loginfo('Before %f'%back_alpha)
            next_back_alpha_commandCode = alpha2commandCode(back_alpha) 
            rospy.loginfo('After %f'%next_back_alpha_commandCode)
            rospy.loginfo("A-::command is %f, %f, next_command is %f, %f  on theta %f, alpha %f, back_alpha %f"%(alpha_commandCode, back_alpha_commandCode, next_alpha_commandCode, next_back_alpha_commandCode, theta, alpha, back_alpha)) 
            rospy.loginfo("m is %f, alpha is %f, %f"%(m_right,alpha, back_alpha)) 
                    #if last_point != 0:
                    #    rospy.loginfo('whats %f'%last_point)
                    #    next_pt = path_lst.pop()
        # slowly update the flipper
            if abs(next_alpha_commandCode - alpha_commandCode) < DELTA_ANGLE_CODE:
                alpha_commandCode = next_alpha_commandCode
            else:
                if next_alpha_commandCode > alpha_commandCode:
                    alpha_commandCode += DELTA_ANGLE_CODE
                else:
                    alpha_commandCode -= DELTA_ANGLE_CODE
            if abs(next_back_alpha_commandCode - back_alpha_commandCode) < DELTA_ANGLE_CODE: 
                back_alpha_commandCode = next_back_alpha_commandCode
            else:
                if next_back_alpha_commandCode > back_alpha_commandCode:
                    back_alpha_commandCode += DELTA_ANGLE_CODE
                else:
                    back_alpha_commandCode -= DELTA_ANGLE_CODE
            rospy.loginfo("A::command is %f, %f, next_command is %f, %f  on theta %f, alpha %f, back_alpha %f"%(alpha_commandCode, back_alpha_commandCode, next_alpha_commandCode, next_back_alpha_commandCode, theta, alpha, back_alpha))
    publish()
    publish_param()
    '''
    command = Float64MultiArray()
    rospy.loginfo("B::command is %f, %f on theta %f"%(alpha_commandCode, back_alpha_commandCode, theta))
    command.data.append(velocity)
    command.data.append(angular) 
    command.data.append(alpha_commandCode)#(front_l)#note that the range is 0~1 that is assume to be pi/2~-pi/2
    command.data.append(alpha_commandCode)#(front_r) 
    command.data.append(back_alpha_commandCode)#(back_l) 
    command.data.append(back_alpha_commandCode)#(back_r) 


    control_pub.publish(command)                                                                                                                                      
    '''
'''        
def publish():
    global velocity, angular, alpha_commandCode, back_alpha_commandCode, control_pub
    command = Float64MultiArray()
    rospy.loginfo("command is %f, %f"%(alpha_commandCode, back_alpha_commandCode))
    command.data.append(velocity)
    command.data.append(angular) 
    command.data.append(alpha_commandCode)#(front_l)#note that the range is 0~1 that is assume to be pi/2~-pi/2
    command.data.append(alpha_commandCode)#(front_r) 
    command.data.append(back_alpha_commandCode)#(back_l) 
    command.data.append(back_alpha_commandCode)#(back_r) 

    control_pub.publish(command) 
'''

def listener():
    global control_pub, record_param_pub
    
    # In ROS, nodes are uniquely named. If two nodes with the same
    # name are launched, the previous one is kicked off. The
    # anonymous=True flag means that rospy will choose a unique
    # name for our 'listener' node so that multiple listeners can
    # run simultaneously.
    rospy.init_node('listener', anonymous=True)
    
        #rospy.Subscriber('/camera/imu', Imu, callback1)
    #rospy.Subscriber('/dynamixel/XM3', XM2, callback2)
    #rospy.Subsriber('/dynamixel/XM4', XM2, callback3)
    #rospy.Subscriber('/vrpn_client_node/camera1/pose', PoseStamped, callback4)

    ts = ApproximateTimeSynchronizer( [Subscriber('/dynamixel/XM3', XM2),
                                       Subscriber('/dynamixel/XM4', XM2)], queue_size=100, slop = 0.1 )
    ts.registerCallback(callback)

    control_pub = rospy.Publisher('Wheel',Float64MultiArray,queue_size=10)
    record_param_pub = rospy.Publisher('record_inner_param',Float64StampedMultiArray,queue_size=10)
                                       
                                       # spin() simply keeps python from exiting until this node is stopped
    rospy.spin()
if __name__ == '__main__':
    #--- navigation and trigger, but now only fix and trigger
    global path_lst, in_planning

    cspace = CSpace(h=H_START, di=D_START)
    path_lst = cspace.path # list of point (d,a,alpha,back_alpha,theta,t)
    path_lst = [cspace.tuple3_to_tuple6(p) for p in path_lst]         
    path_lst = path_lst[::-1]

    in_planning = True
    #---
    try:
        listener()
    except Exception as e:
        print(e)
    finally:
        print('hahahahahahahaha, no problem, I should be on the stair!')
