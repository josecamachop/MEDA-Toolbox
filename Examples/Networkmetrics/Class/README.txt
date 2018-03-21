EDA example with PCA using the MEDA Toolbox. 

Coded by: José Camacho Páez
Last modification of this document: 05/Apr/16

References:

 Data set and Analisys: 

   J. Camacho. Networkmetrics: Multivariate Visual Analytics for Networking Data. 
	Technical Report. 2014.


This folder contains the data and MATLAB scripts for an EDA example based on PCA. 
The data represent network traffic collected through a practical lesson of the Network
Management course of the Telecommunication Engieneering degree at the University of 
Granada. During this session, the students were dealing with the configuration of
polling and traps generation with the Simple Network Management Protocol (SNMP), in the 
Cisco routers and switches of the laboratory.

The Networking Laboratory consists of 24 user work stations (noted Px) connected to the 
different networks configured in the laboratory. These 24 work stations are arranged in 
cells of 4 stations each. The six cells in the laboratory, numbered from 1 to 6, are able 
to work autonomously from the rest of the cells. Each cell is composed by 4 work stations 
(Px/1, Px/2, Px/3 and Px/4, where x stands for the cell number) with three network 
interfaces each one, connected to three different networks containing a variety of network 
devices, including three Cisco 1841 routers (Rx-A, Rx-B and Rx-C) and three Catalyst 2950 
switches (SWx-A, SWx-B and SWx-C), among others (ATM, Frame Relay and X.25 devices, PBX,
etc.).

During part of the laboratory session, SNMP information was acquired from the switches in 
one of the cells (cell 4) with a one minute sampling interval. The data, acquired with the
snmpwalk command of the Net-SNMP distribution, were the input and output traffic octets in 
the 14 interfaces of the three switches. During the session there were two students working 
in the cell: student 1, located in P4/1, was configuring the SW4-A by means of a Telnet 
connection, while student 2, located in P4/2, was configuring R4-A, also with a Telnet 
connection. In addition, during some time intervals in the session, a Neptune attack (SYN
flooding) to SW4-C was performed by the teacher, located in P4/4. The experiment was 
developed so that the students work in the laboratory session was not affected.

The objective of this experiment is to illustrate that the MEDA framework can be used to 
detect anomalies and to identify traffic sources, in particular to perform a forensic 
analysis with the goal of detecting the source of the attack. The experiment is intentionally 
simplistic to simplify the interpretation of EDA for the untrained reader. In this example,
the interpretation of the EDA tools will be performed taking the topology of the network (see
MapasRedLab3.4(A4).pdf) in mind.


Items in the folder:

- data: Original data and software.

- data_proc.mat: Data prepared for analysis. The data set corresponds to 101 observations
	of 49 variables. The list of variables is: 

	1	A.ifIn1		18 	B.ifIn9		35	C.ifIn9
	2	A.ifIn2		19 	B.ifIn10	36	C.ifIn10
	3	A.ifIn8		20 	B.ifIn4		37	C.ifIn12
	4	A.ifIn14	21 	B.ifOut1	38	C.ifIn14
	5	A.ifOut1	22 	B.ifOut2	39	C.ifOut1
	6	A.ifOut2	23 	B.ifOut3	40	C.ifOut2
	7	A.ifOut3	24 	B.ifOut4	41	C.ifOut3
	8	A.ifOut4	25 	B.ifOut8	42	C.ifOut4
	9	A.ifOut8	26 	B.ifOut9	43	C.ifOut5
	10	A.ifOut9	27 	B.ifOut10	44	C.ifOut7
	11	A.ifOut10	28 	B.ifOut11	45	C.ifOut9
	12	A.ifOut11	29 	B.ifOut12	46	C.ifOut10
	13	A.ifOut12	30 	B.ifOut14	47	C.ifOut11
	14	A.ifOut14	31	C.ifIn1		48	C.ifOut12
	15	B.ifIn1		32	C.ifIn2		49	C.ifOut14
	16	B.ifIn2		33	C.ifIn3		
	17	B.ifIn8		34	C.ifIn5		

	where the name of each variable follows the format <switch>.if[In|Out]<N>. These variables
	represent the traffic load (number of octets sent (Out) and received (IN) in a minute in 
	each port (N) of three switches that connect the class network. The remaining ports of the 
	switches did not present any traffic load and where removed form the data.

	Observations have been split into calibration (48 observations) and test (53 obs.)

	The file includes three MATLAB variables:

	- cal: matrix 48x49 with the calibration observations. 

	- test: matrix 53x49 with the test observations.

	- lab: labels of the variables.


- run.m: Script to perform the EDA using the MEDA Toolbox commands. To run the script, the 
	current directory should be the one where the script is stored.

- MapasRedLab3.4(A4).pdf: Diagram with the network topology (in Spanish).
 
- html: web version of the example