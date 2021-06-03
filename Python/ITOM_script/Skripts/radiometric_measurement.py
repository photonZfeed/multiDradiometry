# -*- coding: utf-8 -*-
# =============================================================================
# Created By  : Maximilian Sender
# Copyright   : Copyright 2021, Institute of Chemical Engineering, Prof. Dr. Dirk Ziegenbalg, Ulm University'
# License     : GNU LGPL
# =============================================================================
"""The module radiometric_measurement contains the class ScanSurface, which can be run in an ITOM environment to perform
measurements with a radiometric scanning device developed at the Institute of Chemical Engineering, University Ulm

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np
import pandas as pd
import time
import datetime
import queue
import math
from tinkerforge.ip_connection import IPConnection
from tinkerforge.brick_stepper import BrickStepper
from tinkerforge.bricklet_io16 import BrickletIO16
from tinkerforge.bricklet_industrial_digital_out_4_v2 import BrickletIndustrialDigitalOut4V2
from tinkerforge.brick_silent_stepper import BrickSilentStepper
import matplotlib.pyplot as plt
from Skripts.supplementary.koradserial import KoradSerial  # in case of DLL import error, DLL of the original python3 folder in the itom
# installation should be copied in the working python folder (Anaconda)
import os
from Skripts.supplementary.colour_system import cs_hdtv as cs
from itom import dataObject
from itom import dataIO
from pathlib import Path
from collections import deque
import yagmail


class ScanSurface:
    """
    The class ScanSurface controlling the measurement device is connected to the Spectrometer, the power supply and the
    tinkerforge bricks, which control the stepper motors and end-stop switches. Must be run in an ITOM python
    environment.

    For detailed description see also the documentation of __init__.

        -------------------
        Attributes constructed by __init___:

            reactor (str):
                Name of the reactor (experiment).
            comment (str):
                A comment written to the logfile.
            sendmail (bool):
                Switches on a status update via mail.
            measurements (int):
                Size of the measurement.
            power_supply (bool):
                Switch for power supply (PS) use.
            averages (int):
                Number of averages of spectrometer acquisition.
            x_start (number):
                Scan start in x direction.
            y_start (number):
                Scan start in y direction.
            voltage (number):
                Voltage of PS.
            current (number):
                Current of PS.
            log_output (bool):
                Switch of log file output.
            slicelist:
                List defining the scans z-position.
            int_time_list:
                List defining the integration time of the spectrometer at the respective z-position.
            sendmail (bool, default is False):
                If True, the script will send a status mail after the 50th measurement of a respective z-slice.
                Attributes should be changed in the code.

        -------------------
        Attributes for status and debugging:
            is_homed (bool):
                Indicates over all homing status.
            bool_home_x (bool):
                Indicates x homing status.
            bool_home_y (bool):
                Indicates y homing status.
            z_is_homed (bool):
                Indicates z homing status.
            status (str):
                Indicated status 'idle', 'measuring' or 'homing'.
            path (str):
                The filepath for result export. Normally consisting of '02Results/' the time as "%y%m%d_%H%M%S" and the
                string of the attribute reactor.

        -------------------
        Attributes for the status mail function:
            yag = yagmail.SMTP('fill your email here', 'fill your password here or use the keychain functionality of yagmail')
            mailto (str):
                The receivers e-mail address.
            mailsubj (str):
                Automatically generated subject.
            mailcont (str):
                Automatically generated content.
            approx_string (str):
                Automatically generated string of the approximate ending time of the scan.

        -------------------
        Attributes for Tinkerforge and Steppermotor setup:

            HOST (str):
                "localhost"
            PORT (int):
                4223
            UIDSTEP1 (str):
                UID of Stepper Brick 1
            UIDSTEP2 (str):
                UID of Stepper Brick 2
                Since a core xy-setup is used, there is no such a x or y stepper.
            UIDSTEP3 (str):
                UID of Stepper Brick 3 (z-direction).
            UIDIO (str):
                UID of IO-16 Bricklet.
            UIDIDO4 (str):
                UID of IO-4 Bricklet.
            max_vel (int):
                Velocity of stepper 1 and 2.
            ramping (int):
                Ramping of stepper 1 und 2.
            z_max_vel (int):
                Velocity of z-steppers.
            z_ramping (int):
                Velocity of z-steppers.
            z_ramping (int):
                Ramping of z-steppers.

        -------------------
        Attributes used internally for the homing process:
            goto_queue = queue.Queue()
            spec_queue = queue.Queue()
            home_queue = queue.Queue()
            z_goto_queue = queue.Queue()
            s1_queue = queue.Queue()
            s2_queue = queue.Queue()
            measure_queue = queue.Queue()
            home_queue.put(None)
            home_queue.put(None)
            home_queue.put(None)
        
        -------------------
        Methods: See the respective function's docstring.
            connect()
            start_measurement()
            disconnect()

        Debugging methods for idle mode: (idle mode = connected and homed, no measurement running, mode=='idle')
            get_position()
            get_z_position()
            goto_position(position)
            tune_inttime()
            manual_homing()
            reset_queue()
    """
    def __init__(self,
                 reactor,
                 comment,
                 measurements,
                 power_supply,
                 averages=2,
                 x_start=0.0,
                 y_start=0.0,
                 voltage=5.,
                 current=0.7,
                 log_output=True,
                 linspace=False,  # for setting up the z spacing
                 slicelist=False,
                 int_time_linspace=False,
                 int_time_list=False,
                 sendmail=False
                 ):
        """
        Constructs all the necessary attributes for a ScanSurface instance.
        ------------------
        Parameters:

            reactor (str):
                The name of the experiment, will be used for output folder and files.

            comment (str):
                A comment for additional information, written to the logfile.

            measurements (int):
                The size of measurements of a quadratic canvas in x/y-direction. The points are evenly spaced by 3.9 mm.
                The maximum total size is 153 (59.67 cm).

            power_supply (bool):
                If True, the script controllers the power supply with the given parameters for voltage and current. If
                False, a dark spectrum for the chosen integration time must be provided in the format
                "dark_'integration time'_av'averages'.csv", whereas 'integration time'
                and 'averages' must be replaced with the respective numbers. The provided python script Dunkelspec.py
                can help with that.

            voltage (number, optional, default = 5.0):
                The voltage set for the power supply in Volt.

            current (number, optional, default = 0.7):
                The current set for the power supply in Ampere.

            slicelist (list):
                Sets up a list for all z-positions, in cm, to be measured: [pos1, pos2, pos3...].

            int_time_list (list):
                 Sets up a list of the integration time for all z-positions, in seconds: [time1, time2, time3...].
                 Must be of the same length as parameter slicelist.

            averages (int, optional, default = 2):
                The number of spectrometer measurements which result in one averaged measurement sent to python. A new
                transfer curve is needed for different average settings.

            x_start (number, optional, default = 0.0):
                The x-position of the first measured pixel in cm.

            y_start (number, optional, default = 0.0):
                The y-position of the first measured pixel in cm.

            log_output (boolean, optional, defalut = True):
                Saves a logfile, this is mandatory for further evaluation.

            linspace (tuple, optional, default = False):
                Sets up a slicelist in the manner of the numpy linspace function (start, stop, num). Linspace
                definitions can't be defined together with a list definition (slicelist or int_time_list)!

            int_time_linspace (tuple, optional, default = False):
                Sets up an integration time list in the manner of the numpy linspace function (start, stop, num). Must
                be used together with parameter linspace.
        """
        self.reactor = reactor
        self.comment = comment
        self.sendmail = sendmail
        self.measurements = measurements
        self.power_supply = power_supply
        self.averages = averages
        self.x_start = x_start
        self.y_start = y_start
        self.voltage = voltage
        self.current = current
        self.log_output = log_output
        self.linspace = linspace
        self.slicelist = slicelist
        self.int_time_linspace = int_time_linspace
        self.int_time_list = int_time_list

        if self.linspace is not False and self.slicelist is False:
            if self.int_time_linspace is not False and self.int_time_list is False:
                if len(self.linspace) is 3 and len(self.int_time_linspace) is 3:
                    if self.linspace[2] == self. int_time_linspace[2]:
                        helplist = []
                        self.z_queue = deque(np.linspace(self.linspace[0], self.linspace[1], self.linspace[2]))
                        d1, d2 = self.int_time_linspace[0]**2, self.int_time_linspace[1]**2
                        for i in np.linspace(d1, d2, self.int_time_linspace[2]):
                            helplist.append(round(math.sqrt(i), 0))
                        self.int_time_queue = deque(helplist)
                    else:
                        self.write_out('linspace length not matching')
                        exit()
                else:
                    self.write_out('linspace and int_time_linspace have wrong format')
                    exit()
        elif self.slicelist is not False and self.linspace is False:
            if self.int_time_list is not False and self.int_time_linspace is False:
                if len(self.slicelist) == len(self.int_time_list):
                    self.z_queue = deque(self.slicelist)
                    self.int_time_queue = deque(self.int_time_list)
                else:
                    self.write_out('slicelist and int_time_list are not matching')
                    exit()
        elif self.linspace is not False and self.slicelist is not False:
            print('linspace and slicelist is set, I am using slicelist')
            if len(self.slicelist) == len(self.int_time_list):
                self.z_queue = deque(self.slicelist)
                self.int_time_queue = deque(self.int_time_list)
            else:
                self.write_out('slicelist and int_time_list are not matching')
                exit()

        self.int_time = self.int_time_queue[0]

        self.HOST = "localhost"
        self.PORT = 4223
        self.UIDSTEP1 = "5Wr8fQ"  # UID of Stepper Brick
        self.UIDSTEP2 = "6qY1ds"  # UID of Stepper Brick
        self.UIDSTEP3 = "6wTYWu"  # UID of Stepper Brick
        self.UIDIO = "vYA"  # UID of IO-16 Bricklet
        self.UIDIDO4 = "HjE"
        self.max_vel = 10000
        self.ramping = 10000

        self.z_max_vel = 5000
        self.z_ramping = 10000
        self.z_ramping = 10000

        # homing values:
        self.goto_queue = queue.Queue()
        self.spec_queue = queue.Queue()
        self.home_queue = queue.Queue()
        self.z_goto_queue = queue.Queue()
        self.s1_queue = queue.Queue()
        self.s2_queue = queue.Queue()
        self.measure_queue = queue.Queue()
        self.home_queue.put(None)
        self.home_queue.put(None)
        self.home_queue.put(None)

        self.is_homed = False
        self.bool_home_x = False
        self.bool_home_y = False

        self.z_is_homed = False
        self.bool_home_z = False

        self.should_home = False

        self.status = 'idle'

        self.home_counter = 0
        self.home_counter2 = 0

        self.data = dataObject([1, 2048], 'float32')  # measured intensities

        self.df_for_exp = pd.DataFrame({})
        self.meas_number = False

        self.path = Path(str('02Results/' + time.strftime("%y%m%d_%H%M%S")+'_' + self.reactor))
        self.old_path = self.path
        self.write_out_buffer = []

        self.yag = yagmail.SMTP('fill_your_email@he.re', 'fill your password here or use the keychain functionality of yagmail')
        self.mailto = 'fill the receiving mail address here'
        self.mailsubj = ''
        self.mailcont = ''
        self.approx_string = ''
        self.transfer_object = dataObject([1, 859], 'float32')
        self.dark_data = dataObject([1, 2048], 'float32')  # measured intensities
        self.dark_spec = np.array([])
        self.yaxis = np.array([])
        self.corr_yaxis = np.array([])
        self.xaxis = np.array([])
        self.z_pos_helper = float

        self.integrals = np.array([])
        self.colorarray = np.array([])
        self.x_values = np.array([])
        self.y_values = np.array([])
        self.p_integrals = np.array([])

        self.sat_count = 0
        
    def adapt_lists(self, intlist=False, slicelist=False, couple_z_int=False):
        """Helper method for tune_inttime method."""
        if slicelist is False:
            pass
        elif type(slicelist) is list:
            self.slicelist = slicelist
        else:
            print('slicelist is not a list')
        if type(couple_z_int) is tuple:
            a = couple_z_int[1]/couple_z_int[0]**2
            self.int_time_list = []
            for i in self.slicelist:
                self.int_time_list.append(round((i**2)*a, 3))
        elif couple_z_int is False and type(intlist) is list:
            self.int_time_list = intlist
        print('slicelist:')
        print(self.slicelist)
        print('int time list:')
        print(self.int_time_list)
        print('writing them as queues')
        self.int_time_queue = deque(self.int_time_list)
        self.z_queue = deque(self.slicelist)

    def send_mail(self, key):
        """Sends an e-mail with info about the started measurement and the estimated endtime, called automatically after
        50 measurements"""
        self.mailsubj = 'BLACKBOX STATUS: %s' % key
        self.mailcont = 'Hi there,\n'
        if key == 'ERROR':
            self.mailcont += 'An error occurred. See the attached Log:\n\n'
            for string in self.write_out_buffer:
                if type(string) is str:
                    self.mailcont += (string + '\n')
                elif type(string) is list:
                    for item in string:
                        self.mailcont += (str(item) + ' ')
                    self.mailcont += '\n'
        elif key == 'start':
            self.mailsubj += 'reactor: %s, z position: %.2f' % (self.reactor, self.z_pos_helper)
            self.mailcont += 'New measurement started.\n\n'
            self.mailcont += str("Starting time: " + self.starting_time.strftime('%A, %d. %b %H:%M:%S')+'\n\n')
            self.mailcont += self.approx_string
            self.mailcont += '\n\nThis is the Logfile:\n-------------------------\n\n'
            for string in self.write_out_buffer:
                if type(string) is str:
                    self.mailcont += (string + '\n')
                elif type(string) is list:
                    for item in string:
                        self.mailcont += (str(item) + ' ')
                    self.mailcont += '\n'
        self.yag.send(self.mailto, self.mailsubj, self.mailcont)

    def write_out(self, string):
        """Helper method for writing to the logfile."""
        if self.log_output is False:
            print(string)
        elif self.log_output is True:
            print(string)
            self.write_out_buffer.append(string)

    def change_status(self, string):
        """Helper method for changing the device's status."""
        if string is 'idle' or 'measuring' or 'homing':
            self.write_out('Status changed from ' + self.status + ' to ' + string)
            self.status = string
        elif type(string) is not str:
            print('Wrong type for change_status: ', str(type(string)))
        else:
            print('Wrong keyword for change_status:', string)

    def darkspec(self):
        """Helper method for performing dark measurement or loading a dark spectrum if power _supply is False"""
        print("making dark")
        if self.power_supply is True:
            KoradSerial.device.output.off()
            self.spec.acquire()
            time.sleep(self.int_time*self.averages)
            self.spec.copyVal(self.dark_data)
            self.dark_spec = np.array(self.dark_data[0,203:1062])[0]
            #self.dark_spec = (list(self.dark_data[0, :]))
            #del self.dark_spec[1062:2048]  # adapt measured intensities to calibrated wavelenght range
            #del self.dark_spec[0:203]
            self.df_for_exp['dark'] = self.dark_spec
            print('continue in 5s!!')
            time.sleep(5)
            time.sleep(1)
            KoradSerial.device.output.on()
        else:
            try:
                self.dark_spec = np.genfromtxt(Path(str("Skripts/supplementary/DarkSpectra/dark_"+str(self.int_time).replace('.', '_')+'_av'+str(self.averages)+'.csv'))).T
                self.df_for_exp['dark'] = self.dark_spec
            except IOError:
                raise IOError("Darkfile not found. Check if you created one.")
        print("dark_done")

    def transfer_curve(self):
        """Helper method to load the transfer (calibration) curve."""
        factor_f = (0.03551 / self.int_time)  # integration time of calibration[s]/ integration time measurement [s],
        self.transfer_curve = np.genfromtxt(Path("Skripts/supplementary/TransferCurves/Transfercurve.csv")).T*factor_f
        self.df_for_exp['transfer_curve'] = self.transfer_curve

    def transfer(self, item):
        """Helper method to calculate absolute measures."""
        y_pos = item[1] - self.y_start
        x_pos = item[0] - self.x_start
        if all(i < 63200 for i in self.corr_yaxis) is True:
            pass
        else:
            self.write_out(str('saturation detected at x_pos,y_pos' + str(x_pos) + '/' + str(y_pos)))
            self.sat_count += 1

        transfer_data = self.corr_yaxis * self.transfer_curve
        integral = float(np.trapz(transfer_data, x=self.xaxis))
        p_integral = integral * 0.119  # area of cosine corrector
        integralcc = integral * 0.1521  # scale to a 0.39x0.39 square

        xystr = str(round(x_pos - 0.585, 2)) + '_' + str(round(y_pos - 0.585, 2))
        self.df_for_exp[xystr] = list(transfer_data)
        
        counter = int(
            round(((round(y_pos - 0.585, 2) / 0.39) * self.measurements + (round(x_pos - 0.585, 2) / 0.39)), 0))
        # print((round(y_pos-0.585,2)/0.39)*self.measurements + (round(x_pos-0.585,2)/0.39)+1)
        print('Pixel no: %i  |  %.0f'%(counter,counter/self.measurements**2*100),'%', end='\r')
        self.integrals[counter] = integralcc
        self.p_integrals[counter] = p_integral
        self.x_values[counter] = round(x_pos, 3)
        self.y_values[counter] = round(y_pos, 3)

        # calculate color, self.colorarray was set up before
        lam = np.arange(380., 781., 5)  # lambda table for spec_to_xyz
        self.colorarray[counter] = cs.spec_to_rgb(np.interp(lam, self.xaxis, transfer_data))

    def spec_measure(self, x):
        """Helper method, queuing a spectrometer measurement."""
        self.spec_queue.put([x])

    def measure(self, item):
        """Helper method calling a spectrometer measurement, for end time estimation and counting measurements."""
        item2 = self.spec_queue.get()
        if item2 is None:
            self.write_out("item is None, no measurement")
        else:
            self.spectrometer(item)
            if self.meas_number is not False:
                self.meas_number += 1
            else:
                self.meas_number = 0
            if self.meas_number is 50:
                self.approx_string = str("Approx end time: " + (self.starting_time3 + datetime.timedelta(
                    seconds=(((time.time() - self.starting_time2) / 50) * (self.measurements ** 2)))).strftime('%A, %d. %b %H:%M:%S'))
                self.write_out(self.approx_string)
                if self.sendmail is True:
                    self.send_mail('start')
            else:
                pass

    def spectrometer(self, item):
        """Helper method to perform the spectrometer measurement."""
        self.spec.acquire()
        time.sleep(self.int_time*self.averages)
        self.spec.copyVal(self.data)
        self.yaxis = np.array(self.data[0,203:1062])[0]
        self.corr_yaxis = self.yaxis - self.dark_spec
        if self.meas_number is False:
            self.xaxis = np.array(self.spec.getParam("lambda_table"))[203:1062]
            self.df_for_exp['wavelength'] = self.xaxis
        self.transfer(item)
        self.measure_queue.task_done()

    def get_position(self):
        """Gives the current position of the slid as list with size 2, in cm"""
        x_pos = (0.5 * (self.stepper1.get_current_position() + self.stepper2.get_current_position())) / 333.33
        y_pos = (0.5 * (self.stepper1.get_current_position() - self.stepper2.get_current_position())) / 333.33
        return [x_pos, y_pos]

    def get_z_position(self):
        """Gives the current position of the z_axis in cm"""
        # 1600 steps per rev
        # 3mm per rev
        # 533 1/3 steps per rev
        return self.stepper3.get_current_position()*(0.3/1600)

    def put_goto_position(self, x_pos, y_pos):
        """Helper method to queue spectrometer position for a measurement."""
        self.goto_queue.put([x_pos, y_pos])

    def goto_z_position(self, z_pos, measure=False):
        """drives the z slid in the given position in cm, if measure is True, the pos reached of the z stepper will start
        the measurement"""
        self.z_pos_helper = z_pos
        if z_pos < 0 or z_pos > 55:
            self.write_out('z position outside boundaries')
            pass
        else:
            if measure is True:
                self.path = self.old_path / str(str(z_pos) + '_mm')
            self.z_goto_queue.put(None)
            self.stepper3.set_target_position(round(-1*z_pos*(6400/0.3)))  # NEU 200*32 und target pos

    def z_reached(self, x):
        """Helper method, called when z-endstop reached."""
        time.sleep(1)
        self.z_goto_queue.task_done()
        print('z_arrived')

    def worker(self, x):
        """Helping method, sends the slid to the position, then waits for it to arrive
        and starts the measurement, then is waiting for the measurement to finish."""
        if self.status is 'measuring':
            working = True
            while working is True:
                item = self.goto_queue.get()
                if item is None:
                    self.write_out('item is None')
                    self.goto_queue.task_done()
                elif item is 'end':
                    working = False
                    self.goto_position_inqueue([0.2, 0.2])
                    self.s1_queue.join()
                    self.s2_queue.join()
                    self.goto_queue.task_done()
                elif item[0] == 0 and item[1] == 0:
                    self.goto_position_inqueue(item)
                    self.goto_queue.task_done()
                else:
                    self.goto_position_inqueue(item)
                    self.s1_queue.join()
                    self.s2_queue.join()
                    self.measure_queue.put('None')
                    self.measure(item)
                    self.measure_queue.join()
                    self.goto_queue.task_done()
        elif self.status is 'homing':
            pass
        elif self.status is 'idle':
            pass

    def goto_position_inqueue(self, item):
        """Helper method. Sends the slid to the given position, in cm. Only to be called by the worker method."""
        self.s1_queue.put(None)
        self.s2_queue.put(None)
        x_pos = item[0]
        y_pos = item[1]
        step_pos_1 = 333.33 * (x_pos + y_pos)
        step_pos_2 = 333.33 * (x_pos - y_pos)
        self.stepper1.set_steps(int(step_pos_1 - self.stepper1.get_current_position()))
        self.stepper2.set_steps(int(step_pos_2 - self.stepper2.get_current_position()))
        if round((x_pos - self.x_start), 2) == 0.39 and round((y_pos - self.y_start), 2) == 0.39:
            pass
        else:
            pass

    def goto_position(self, item):
        """Drives the spectrometer to the given position. Works only if mode is 'idle'.
            Parameters:
                item: (list), [x-position, y-position] in cm."""
        if self.status is 'idle':
            self.stepper1.set_max_velocity(2000)
            self.stepper2.set_max_velocity(2000)
            x_pos = item[0]
            y_pos = item[1]
            step_pos_1 = 333.33 * (x_pos + y_pos)
            step_pos_2 = 333.33 * (x_pos - y_pos)
            self.stepper1.set_steps(int(step_pos_1 - self.stepper1.get_current_position()))
            self.stepper2.set_steps(int(step_pos_2 - self.stepper2.get_current_position()))
            time.sleep(0.1)
        else:
            print('not in idle mode')

    def cb_interrupt(self, port, interrupt_mask, value_mask):
        """Helper method, called when an endstop is triggered."""
        print("something triggered")
        print("value mask: ", str(value_mask))
        if self.status is 'homing':
            if value_mask == 248:  # 11111000
                if self.bool_home_z is True:
                    self.stepper3.full_brake()
                    self.stepper3.disable()
                    self.stepper3.enable()
                    self.stepper3.set_current_position(0)
                    self.bool_home_z = False
                    self.z_is_homed = True
                    self.write_out('z homed')
                    self.home_queue.task_done()
                    print('q_size: ' + str(self.home_queue.qsize()))
                else:
                    self.stepper3.disable()
                    self.write_out('Z triggered during measurement, stopping')
            if value_mask == 252 and self.bool_home_x is True:  # x: 254, y: 253, both: 252 mask = 11111100
                self.stepper1.full_brake()
                self.stepper2.full_brake()
                self.write_out(
                    str('Endstop x triggered, motor stopped, self.current pos:' + str(self.stepper1.get_current_position())))
                if self.bool_home_x is True and self.home_counter2 == 0:
                    self.write_out('home mode, first contact, retracting')
                    self.stepper1.set_max_velocity(2000)
                    self.stepper2.set_max_velocity(2000)
                    self.stepper1.set_steps(160)
                    self.stepper2.set_steps(160)
                elif self.bool_home_x is True and self.home_counter2 >= 1:
                    self.bool_home_x = False
                    self.write_out('driving home')
                    time.sleep(1)
                    self.stepper1.set_current_position(0)
                    self.stepper2.set_current_position(0)
                    self.write_out(str('self.current pos1:' + str(self.stepper1.get_current_position())))
                    self.write_out(str('self.current pos2:' + str(self.stepper2.get_current_position())))
                    self.write_out('--------------------homing process x finished--------------------')
                    self.home_queue.task_done()
                    time.sleep(1)
                    self.should_home = True
                else:
                    pass

            if value_mask == 253 and self.bool_home_y is True:  # x: 254, y: 253, both: 252 mask = 11111101
                self.stepper1.full_brake()
                self.stepper2.full_brake()
                self.write_out(
                    str('Endstop y triggered, motor stopped, self.current pos:' + str(self.stepper1.get_current_position())))
                if self.bool_home_y is True and self.home_counter == 0:
                    self.write_out('home mode, first contact, retracting')
                    self.stepper1.set_max_velocity(2000)
                    self.stepper2.set_max_velocity(2000)
                    self.stepper1.set_steps(160)
                    self.stepper2.set_steps(-160)
                if self.bool_home_y is True and self.home_counter > 0:
                    self.bool_home_y = False
                    self.write_out('driving home')
                    time.sleep(1)
                    self.stepper1.set_current_position(0)
                    self.stepper2.set_current_position(0)
                    self.write_out('--------------------homing process y finished--------------------')
                    self.home_queue.task_done()
                    time.sleep(1)
                    self.should_home = True
                    # self.stepper2.disable()
                    # self.ipcon.disconnect
            elif self.bool_home_x is True and value_mask == 253:  # mask = 11111101
                self.write_out(str('X Endstop open while homing (' + str(self.home_counter2) + ')'))
                time.sleep(1)
                self.home_counter2 += 1
                self.stepper1.set_max_velocity(1000)
                self.stepper2.set_max_velocity(1000)
                self.stepper1.set_steps(-160)
                self.stepper2.set_steps(-160)

            if self.bool_home_y is True and value_mask == 255:  # mask = 11111111
                self.write_out(str('Y Endstop open while homing (' + str(self.home_counter) + ')'))
                time.sleep(1)
                self.home_counter += 1
                self.stepper1.set_max_velocity(1000)
                self.stepper2.set_max_velocity(1000)
                self.stepper1.set_steps(-160)
                self.stepper2.set_steps(160)
            elif self.bool_home_y is False and self.bool_home_z is True and value_mask == 255:
                self.bool_home_z = False  # z-endstop opened for x and y homing, z will home after that
                self.write_out('Z Endstop was triggered, opening, then proceeding to home Y')
                self.stepper3.full_brake()
                self.should_home = True
            elif self.bool_home_z is False and self.bool_home_y is False and self.bool_home_x is False and value_mask == 255:  # 11111111
                self.write_out(str('Z Endstop opened'))

            if self.should_home is True:
                self.should_home = False
                self.home()
        elif self.status is not 'homing':
            self.write_out('Endstop hit while Status' + self.status)

    def reset_home_queue(self):
        """Resets the home-queue in order to home an already homed instance again."""
        if self.home_queue.qsize is 3:
            raise Exception('home_queue is already 3')
        else:
            for i in range(0, self.home_queue.qsize()):
                self.home_queue.task_done()
            self.home_queue.put(None)
            self.home_queue.put(None)
            self.home_queue.put(None)
            self.write_out('home_queue is set to 3')

    def home(self):
        """Internal function! Gets called by start_measurement"""
        if self.home_queue.qsize() <= 3 > 0:
            pass
        elif self.home_queue.qsize() == 0:
            raise Exception('Already homed. Do you want to home again? Then reset the home queue by reset_home_queue().')

        print('home function sees: ', str(self.io16.get_port("a")))
        self.stepper1.set_max_velocity(2000)
        self.stepper2.set_max_velocity(2000)
        self.stepper1.enable()
        self.stepper2.enable()
        if self.io16.get_port("a") is 255:  # 11111111
            self.write_out('--------------------homing y process started--------------------')
            self.bool_home_y = True
            self.stepper1.set_steps(-100000)
            self.stepper2.set_steps(100000)
        elif self.io16.get_port("a") is 253:  # 11111101
            self.write_out('--------------------homing x process started--------------------')
            self.bool_home_x = True
            self.stepper1.set_steps(-100000)
            self.stepper2.set_steps(-100000)
        elif self.io16.get_port("a") is 252:  # 11111100
            self.write_out('--------------------homing z process started--------------------')
            self.bool_home_z = True
            self.stepper3.drive_forward()
        elif self.io16.get_port("a") is 248:  # 11111000
            self.write_out("homing finished")
        elif self.io16.get_port("a") is 251:  # 11111011
            print('only z triggered but not homed, driving backwards')
            self.bool_home_z = True
            self.stepper3.set_steps(-21000)
        else:
            print('should home, but nothing applies')

    def signaltonoise(self, a, axis=0, ddof=0):
        """Helper method, estimating the signals quality."""
        a = np.asanyarray(a)
        m = a.mean(axis)
        print('mean: ', str(m))
        sd = a.std(axis=axis, ddof=ddof)
        print('std: ', str(sd))
        return np.where(sd == 0, 0, m/sd)

    def tune_inttime(self):
        """Tunes the inttime to an appropriate time, depending on the measured intensity."""
        self.write_out('Tuning integration time.')
        hit = False
        run = 1
        itime = self.int_time
        while hit is False and 0.003 < itime < 3:
            self.spec.setParam("integration_time", itime)
            self.spec.acquire()
            time.sleep(itime*self.averages)
            self.spec.copyVal(self.data)

            values = np.array(self.data[0,203:1062])[0]
            
            self.write_out('run number '+str(run))

            self.write_out('integration time: '+str(itime))
            self.write_out('signal to noise ratio: '+str(self.signaltonoise(values)))
            self.write_out('signal maximum: '+str(max(values)))
            self.write_out('-----------------------------------------')
            if max(values) > 63200:
                print('over 63200')
                itime = itime/2
            elif max(values) < 63200:
                print('under 63200')
                if self.signaltonoise(values) < 2:

                    if 50000 < max(values) < 56000:
                        hit = True
                    else:
                        itime = itime*(55000/max(values))
                else:
                    itime = itime*2
            if run == 100:
                hit = True
            run += 1
        if hit is True and run != 101:
            self.int_time = itime
            self.write_out('Found a matching integration time of '+str(itime)+' seconds')
            user = input('should I adapt the int time list according to it? y/N')
            if user == 'y' or 'Y':
                self.write_out('Adapting int time list')
                self.adapt_lists(couple_z_int=(self.z_pos_helper, itime))
        elif run == 101:
            self.write_out('could not find signal after 100 iterations')
        elif 0.003 < itime < 3:
            self.write_out('could not find signal inside integration time boundaries')

    def manual_homing(self):
        """Starts homing procedure without starting the measurement."""
        print('im homing')
        self.change_status('homing')
        self.home()
        self.home_queue.join()
        self.write_out("homing ended")
        self.stepper1.set_max_velocity(self.max_vel)
        self.stepper2.set_max_velocity(self.max_vel)
        self.stepper3.register_callback(self.stepper3.CALLBACK_POSITION_REACHED, lambda x: self.z_reached(x))
        """setting the callbacks of the steppers to a function which ends the job in the queue which was raised by the 
        order to go somewhere"""
        self.stepper1.register_callback(self.stepper1.CALLBACK_POSITION_REACHED, lambda x: self.s1arrived(x))
        self.stepper2.register_callback(self.stepper1.CALLBACK_POSITION_REACHED, lambda x: self.s2arrived(x))
        self.is_homed = True
        self.change_status('idle')

    def s1arrived(self, x):
        self.s1_queue.task_done()

    def s2arrived(self, x):
        self.s2_queue.task_done()

    def start_measurement(self):
        """Starts the measurement. Homing will be performed if the device is not homed."""
        self.int_time = self.int_time_queue.popleft()
        self.spec.setParam("integration_time", self.int_time)
        self.spec.setParam("average", self.averages)
        self.spec.setParam("dark_correction", 2)  # dyndark enabled
        self.meas_number = False
        self.write_out(str("Name: " + self.reactor))
        self.write_out(self.comment)
        self.write_out("Size: " + str(self.measurements))
        self.write_out("x_start: " + str(self.x_start))
        self.write_out("y_start: " + str(self.y_start))
        self.write_out("int_time: " + str(self.int_time))
        self.write_out("averages: " + str(self.averages))
        self.write_out("Power_supply: " + str(self.power_supply))

        if self.power_supply is True:
            self.write_out("voltage: " + str(self.voltage))
            self.write_out("current: " + str(self.current))
        else:
            pass

        if self.power_supply is True:
            """Only switches on power supply, connect is in connect function"""
            self.m1.recall()
            self.channel.voltage = self.voltage
            self.channel.current = self.current
            self.m1.save()
            self.m1.recall()
            KoradSerial.device.output.on()

        # lists for results:
        self.integrals = np.zeros((self.measurements ** 2))
        self.colorarray = np.zeros((self.measurements ** 2, 3))
        self.x_values = [0] * self.measurements ** 2
        self.y_values = [0] * self.measurements ** 2
        self.p_integrals = [0] * self.measurements ** 2

        self.sat_count = 0

        self.darkspec()
        time.sleep(3)
        self.transfer_curve()

        if self.is_homed is False:
            print('im homing')
            self.change_status('homing')
            self.home()
            self.home_queue.join()
            self.write_out("homing ended")
            self.stepper1.set_max_velocity(self.max_vel)
            self.stepper2.set_max_velocity(self.max_vel)
            # self.stepper1.set_current_position(0)
            # self.stepper2.set_current_position(0)
            self.stepper3.register_callback(self.stepper3.CALLBACK_POSITION_REACHED, lambda x: self.z_reached(x))

            """setting the callbacks of the steppers to a function which ends the job in the queue which was raised by the 
            order to go somewhere"""
            self.stepper1.register_callback(self.stepper1.CALLBACK_POSITION_REACHED, lambda x: self.s1arrived(x))
            self.stepper2.register_callback(self.stepper1.CALLBACK_POSITION_REACHED, lambda x: self.s2arrived(x))
            self.is_homed = True
        self.change_status('idle')
        self.z_pos_helper = self.z_queue.popleft()
        self.write_out("z_pos: " + str(self.z_pos_helper))
        self.goto_z_position(self.z_pos_helper, True)
        self.z_goto_queue.join()

        t = 0.195 + self.y_start
        for y_direction in range(0, int(self.measurements / 2)):  # number of self.measurements in y-direction
            r = 0.39 + 0.195 + self.x_start
            t += 0.39
            self.put_goto_position(r, t)
            self.spec_measure(1)
            for x_direction in range(0, (self.measurements - 1)):  # number of self.measurements in x-direction
                r += 0.39
                self.put_goto_position(r, t)
                self.spec_measure(1)
            t += 0.39
            self.put_goto_position(r, t)
            self.spec_measure(1)
            for x_direction in range(0, self.measurements - 1):  # number of self.measurements in x-direction
                r -= 0.39
                self.put_goto_position(r, t)
                self.spec_measure(1)
        if self.measurements % 2 == 1:
            t += 0.39
            self.put_goto_position(r, t)
            self.spec_measure(1)
            for x_direction in range(0,
                                     self.measurements - 1):  # number of self.measurements in y-direction, for uneven size
                r += 0.39
                self.put_goto_position(r, t)
                self.spec_measure(1)
        else:
            pass

        self.starting_time = datetime.datetime.now()
        self.starting_time2 = time.time()
        self.starting_time3 = datetime.datetime.now()
        self.write_out(str("Strating time: " + self.starting_time.strftime('%A, %d. %b %H:%M:%S')))

        self.change_status('measuring')
        self.stepper1.set_max_velocity(self.max_vel)
        self.stepper2.set_max_velocity(self.max_vel)

        """this is the new end method, the worker function gets the keyword end, drives to 0,2 0,2 , does NOT 
        perform a measurement, sets the goto queue and breaks the while loop of worker"""
        self.goto_queue.put('end')
        self.write_out('all put')

        self.worker('start')

        """with join(), we wait for all measurements to finish"""
        self.goto_queue.join()

        self.change_status('idle')

        self.write_out(str("saturations detected: " + str(self.sat_count)))
        self.write_out(str("start: " + self.starting_time.strftime('%A, %d. %b %H:%M:%S') + "   end: " + time.strftime('%A, %d. %b %H:%M:%S')))
        self.write_out(str("This took " + str(datetime.timedelta(seconds=int(time.time() - self.starting_time2))) + " hours"))
        if self.power_supply is True:
            KoradSerial.device.output.off()
        self.export()

    def callback_dummy(self, x):
        pass

    def connect(self):
        """Connects the class instance to the hardware."""
        self.ipcon = IPConnection()  # Create IP connection
        self.stepper1 = BrickStepper(self.UIDSTEP1, self.ipcon)  # Create device object
        self.stepper2 = BrickStepper(self.UIDSTEP2, self.ipcon)  # Create device object

        print("stepper3 connected")
        self.stepper3 = BrickSilentStepper(self.UIDSTEP3, self.ipcon)  # Create device object
        self.io16 = BrickletIO16(self.UIDIO, self.ipcon)  # Create device object

        self.ido4 = BrickletIndustrialDigitalOut4V2(self.UIDIDO4, self.ipcon)

        self.ipcon.connect(self.HOST, self.PORT)  # Connect to brickd
        # Don't use device before self.ipcon is connected
        try:
            self.spec = dataIO("AvantesAvaSpec", 6546, 1639)
        except RuntimeError:
            del self.spec
            self.spec = dataIO("AvantesAvaSpec", 6546, 1639)
        self.spec.startDevice()
        self.spec.startDevice()
        self.spec.setParam("integration_time", self.int_time)
        self.spec.setParam("average", self.averages)

        self.ido4.set_pwm_configuration(0, 120000, 700)

        self.stepper1.disable()
        self.stepper2.disable()
        self.stepper3.disable()
        self.stepper1.enable()
        self.stepper2.enable()
        self.stepper3.enable()

        self.stepper1.set_motor_current(1200)  # mA
        self.stepper1.set_step_mode(8)  # 1/8 step mode
        self.stepper1.set_max_velocity(self.max_vel)  # Velocity: steps/s
        self.stepper1.set_speed_ramping(self.ramping, self.ramping)
        self.stepper1.set_sync_rect("true")
        self.stepper1.set_decay(60000)

        self.stepper2.set_motor_current(1200)  # mA
        self.stepper2.set_step_mode(8)  # 1/8 step mode
        self.stepper2.set_max_velocity(self.max_vel)  # Velocity: steps/s
        self.stepper2.set_speed_ramping(self.ramping, self.ramping)
        self.stepper2.set_sync_rect("true")
        self.stepper2.set_decay(60000)

        self.stepper3.set_motor_current(800)  # 800mA
        self.stepper3.set_step_configuration(self.stepper3.STEP_RESOLUTION_32, True)  # 1/8 steps (interpolated)
        self.stepper3.set_speed_ramping(20000, 65535)
        self.stepper3.set_max_velocity(13000)
        self.stepper3.register_callback(self.stepper3.CALLBACK_POSITION_REACHED, lambda x: self.z_reached(x))

        # Enable interrupt on pin 0 of port a
        self.io16.register_callback(self.io16.CALLBACK_INTERRUPT, self.cb_interrupt)
        self.io16.set_port_configuration("a", 255, "i", True)
        self.io16.set_port_interrupt("a", 7)

        self.stepper1.set_max_velocity(1000)
        self.stepper2.set_max_velocity(1000)
        self.stepper1.set_steps(1600)

        if self.power_supply is True:
            """Only connects to power supply,switch on is in start function"""
            KoradSerial.device = KoradSerial('COM5', True)
            self.write_out("LED connected")
            self.channel = KoradSerial.device.channels[0]
            self.m1 = KoradSerial.device.memories[0]

        time.sleep(3)
        self.stepper1.set_speed_ramping(self.ramping, self.ramping)

    def export(self):
        """Exports the full measurement in the folder "02Results", automatically called."""
        if not os.path.exists(self.path):
            os.makedirs(self.path)  # create new folder containing all results
        else:
            print("Folder already exists! Writing to", self.path, '.')

        with open(self.path / str(self.reactor + '.log'), 'a') as f:
            for string in self.write_out_buffer:
                if type(string) is str:
                    f.write(string + '\n')
                elif type(string) is list:
                    with open(self.path / str(self.reactor + '.log'), 'a') as f:
                        for item in string:
                            f.write(str(item) + ' ')
                        f.write('\n')

        planck_constant = 6.626 * 10 ** (-34)
        speed_of_light = 2.998 * 10 ** 8
        avogadro_number = 6.022 * 10 ** 23

        p_results = list(
            np.array(self.p_integrals) * 10 ** (-9) / planck_constant / speed_of_light / avogadro_number
            * 10 ** (-6))

        plt.ioff()
        self.df_for_exp.to_feather(self.path / str(self.reactor + ".feather"))
        x, y = np.meshgrid(np.linspace(0., 0.39 * (self.measurements+1), self.measurements+1), np.linspace(0., 0.39 * (self.measurements+1), self.measurements+1))
        z = np.reshape(np.array(self.integrals), (self.measurements, self.measurements))
        plt.xlabel("x / cm")
        plt.ylabel("y / cm")
        plt.pcolormesh(x, y, z, antialiased=False)
        plt.colorbar(label="P / \u03BCW")
        plt.title('received power: ' + str(np.round(np.sum(self.integrals)*1E-6, 3)) + ' W')
        plt.axis([x.min(), x.max(), y.min(), y.max()])
        plt.xticks(np.linspace(0., 0.39 * (self.measurements+1), 7))
        plt.yticks(np.linspace(0, 0.39 * (self.measurements+1), 7))
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig(str(self.path / str('Colormap' + self.reactor + '.png')), dpi=300, bbox_inches='tight')
        plt.clf()

        """Plot Color"""
        x, y = np.meshgrid(np.linspace(0., 0.39 * (self.measurements+1), self.measurements+1), np.linspace(0., 0.39 * (self.measurements+1), self.measurements+1))
        z = np.reshape(np.array(self.colorarray), (self.measurements, self.measurements, 3))
        plt.xlabel("x / cm")
        plt.ylabel("y / cm")
        plt.imshow(z, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
        plt.xticks(np.linspace(0., 0.39 * (self.measurements+1), 7))
        plt.yticks(np.linspace(0, 0.39 * (self.measurements+1), 7))
        plt.savefig(str(self.path / str('Farbbild' + self.reactor + '.png')), dpi=300, bbox_inches='tight')
        plt.clf()
        plt.close()
        # if sliced and que not empty start measurement
        if len(self.z_queue) is not 0:
            if self.power_supply is True:
                KoradSerial.device.output.off()
            self.start_measurement()
        else:
            pass

    def plot_fig_rad(self):
        """Plotting"""
        plt.ion()
        x, y = np.meshgrid(np.linspace(0., 0.39 * (self.measurements+1), self.measurements+1), np.linspace(0., 0.39 * (self.measurements+1), self.measurements+1))
        z = np.reshape(np.array(self.integrals), (self.measurements, self.measurements))
        plt.xlabel("x / cm")
        plt.ylabel("y / cm")
        plt.pcolormesh(x, y, z, antialiased=False)
        plt.colorbar(label="P / \u03BCW")
        plt.title('received power: ' + str(np.round(np.sum(self.integrals)*1E-6, 3)) + ' W')
        plt.axis([x.min(), x.max(), y.min(), y.max()])
        plt.xticks(np.linspace(0., 0.39 * (self.measurements+1), 7))
        plt.yticks(np.linspace(0, 0.39 * (self.measurements+1), 7))
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    def plot_fig_color(self):
        """Plotting"""
        plt.ion()
        """Plot Color"""
        """Plotting"""
        self.colorarray = np.clip(self.colorarray, 0, 1)
        x, y = np.meshgrid(np.linspace(0., 0.39 * (self.measurements+1), self.measurements+1), np.linspace(0., 0.39 * (self.measurements+1), self.measurements+1))
        z = np.reshape(np.array(self.colorarray), (self.measurements, self.measurements, 3))
        plt.xlabel("x / cm")
        plt.ylabel("y / cm")
        plt.imshow(z, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
        plt.xticks(np.linspace(0., 0.39 * (self.measurements+1), 7))
        plt.yticks(np.linspace(0, 0.39 * (self.measurements+1), 7))
        plt.show()

    def disconnect(self):
        """Disconnects the script from the hardware."""
        self.stepper1.disable()
        self.stepper2.disable()
        self.stepper3.disable()
        self.ido4.set_pwm_configuration(0, 0, 0)
        self.ido4.set_value([False, False, False, False])
        self.ipcon.disconnect()
        self.spec.stopDevice()
        del self.spec
        if self.power_supply is True:
            KoradSerial.device.output.off()
            KoradSerial.device.close()


__author__ = "Maximilian Sender"
__copyright__ = 'Copyright 2021, Institute of Chemical Engineering, Ulm University'
__license__ = "GNU LGPL"
