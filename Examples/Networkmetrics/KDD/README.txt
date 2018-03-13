EDA example for Big Data using the MEDA Toolbox.

Coded by: José Camacho Páez
Last modification of this document: 05/Apr/16

References:

 Data set and Analisys: 

   J. Camacho. Networkmetrics: Multivariate Visual Analytics for Networking Data. 
	Technical Report. 2014.

   J. Camacho. Visualizing Big data with Compressed Score Plots: Approach and Research 
    	Challenges. Chemometrics and Intelligent Laboratory Systems, vol. 135, pp. 
	110-125, 2014.


the data set was generated from the 1998 DARPA Intrusion Detection evaluation Program, 
prepared and managed by MIT Lincoln Labs. The objective of this program was to survey 
and evaluate research in networking intrusion detection. For that, a large data set 
including a wide variety of intrusions simulated in a military network environment
was provided. The original data set included 4.880.000 observations (connection records). 
The observations belong to 22 different classes, one class for normal traffic and the 
remaining for different types of network attacks. Four main categories of attacks were 
simulated: denial-of-service (DoS), e.g., syn flood; unauthorized access from a remote 
machine, e.g., guessing password; unauthorized access to local superuser (root) privileges, 
e.g., buffer overflow attacks; and surveillance and probing, e.g., port-scan. For 
illustrative purposes, the analysis will be restricted to two types of DoS attack, smurf 
and neptune, and normal traffic. These three classes represent a 99.3% of the total traffic 
in the data set. For each connection, 42 features are computed, including numerical and 
categorical features. To consider categorical features in the EDA, one dummy variable per 
category is included in the data set. The resulting data set contains 4.844.253 observations, 
each one with 122 features.

To analyze Big Data, the data should be stored in a set of files and the list with the 
files name should be passed as a function argument.

The variables in the data are:

- Numerical

'dur', %'duration',
'srcd', %'src_bytes',
'dstb', %'dst_bytes',
'land', %'land',
'wf', %'wrong_fragment',
'urg', %'urgent',
'hot', %'hot',
'nfl', %'num_failed_logins',
'li', %'logged_in',
'nc', %'num_compromised',
'rs', %'root_shell',
'sa', %'su_attempted',
'nr', %'num_root',
'nfc', %'num_file_creations',
'ns', %'num_shells',
'naf', %'num_access_files',
'noc', %'num_outbound_cmds',
'ihl', %'is_host_login',
'igl', %'is_guest_login',
'cnt', %'count',
'sc', %'srv_count',
'sr', %'serror_rate',
'ssr', %'srv_serror_rate',
'rr', %'rerror_rate',
'srr', %'srv_rerror_rate',
'ssr2', %'same_srv_rate',
'dsr', %'diff_srv_rate',
'sdhr', %'srv_diff_host_rate',
'dhc', %'dst_host_count',
'dhsc', %'dst_host_srv_count',
'dhssr', %'dst_host_same_srv_rate',
'dhdsr', %'dst_host_diff_srv_rate',
'dhsspr', %'dst_host_same_src_port_rate',
'dhsdhr', %'dst_host_srv_diff_host_rate',
'dhsr', %'dst_host_serror_rate',
'dhssr', %'dst_host_srv_serror_rate',
'dhrr', %'dst_host_rerror_rate',
'dhsrr', %'dst_host_srv_rerror_rate'

- Binary developed from categorical

'pt0-2', % prottype
'srv0-69', % service
'fg0-10', % flag

Items in the folder:

- Data: Data prepared for analysis.

- kdd.mat: List with the names of the files where the data are stored.
	- list: 489 data files (it takes hours to process)
	- short_list: 10 data files, suitable for the example
	- label: labels of the variables

- run.m: Script to perform the EDA using the MEDA Toolbox commands. To run the script, the 
	current directory should be the one where the script is stored.
 
- html: web version of the example
