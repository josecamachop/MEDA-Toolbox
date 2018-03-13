EDA example with GPCA using the MEDA Toolbox. 

Coded by: José Camacho Páez
Last modification of this document: 25/Jan/17

References:

 Data set and Analisys: 

   J. Camacho, R.A. Rodríguez-Gómez, E. Saccenti, Group-wise Principal Component Analysis 
	for Exploratory Data Analysis. Accepted in Journal of Computational and Graphical 
	Statistics, 2017.


The VAST 2012 2nd mini challenge is a benchmark for visualization in cybersecurity 
(http://www.vacommunity.org/VAST+Challenge+2012)

The goal is to identify cybersecurity issues in the data collected during two days from a
computer network. During those days, a number of non-legitimate programs were found
to be running on several computers, slowing them down. A cyber-forensics operation is
required to discover the root causes for this strange behavior.

Two typical sources of data are collected from the network: firewall and Intrusion
Detection System (IDS) logs. The firewall analyses the incoming and outgoing data traffic
in the network, and records in a log file all connection attempts that are blocked according
to security policies. The IDS employs higher level intelligence to identify cybersecurity
incidents in data traffic. It also stores the results in a log file, though it does not block any
trac connection. Also, typically, it only analyses a sub-set (sample) of the total traffic.

A total of 2350 observations, each one with the information for one minute, are obtained.
For every sampling period of one minute, we defined a set of 112 variables that represent
the information from the two data sources: 69 variables for the firewall log and 43 for
the IDS log. The number of variables was reduced to 95 by discarding those with constant 
value throughout the capture period. 

The variables in the data are:

    'fw_asa66'
    'fw_ipfwhq'
    'fw_empty'
    'fw_syserror'
    'fw_pnetbios3'
    'fw_syscritical'
    'fw_outbound'
    'fw_conbuilt'
    'fw_pshell'
    'fw_ipweb'
    'fw_asa611'
    'fw_asa610'
    'fw_asa639'
    'fw_inbound'
    'fw_pnstd'
    'fw_asa635'
    'fw_asa634'
    'fw_asa636'
    'fw_asa630'
    'fw_asa633'
    'fw_asa671'
    'fw_asa673'
    'fw_asa672'
    'fw_pftp'
    'fw_asa677'
    'fw_sysnotice'
    'fw_opteardown'
    'fw_syswarn'
    'fw_opdeny'
    'fw_opcommand'
    'fw_phttp'
    'fw_pdns'
    'fw_opbuilt'
    'fw_pnetbios2'
    'fw_udp'
    'fw_psmb'
    'fw_phttps'
    'fw_asa21'
    'fw_pnetbios1'
    'fw_psnmp'
    'fw_ptelnet'
    'fw_asa6310'
    'fw_iplog'
    'fw_ipdns'
    'fw_asa33'
    'fw_asa37'
    'fw_asa4'
    'fw_asa5'
    'fw_asa6'
    'fw_pdce'
    'fw_asa3'
    'fw_ipdc'
    'fw_sysinfo'
    'fw_pldap'
    'fw_ipws'
    'fw_opempty'
    'fw_asa517'
    'fw_asa518'
    'fw_condown'
    'fw_ipfwr'
    'fw_pident'
    'fw_ipout'
    'fw_pssh'
    'fw_pkerberos'
    'fw_asa44'
    'fw_tcp'
    'fw_asa41'
    'fw_opdenyacl'
    'ids_privacy'
    'ids_lsnmp'
    'ids_psmb'
    'ids_ipws'
    'ids_lvnc'
    'ids_badtraffic'
    'ids_pdns'
    'ids_ipfwhq'
    'ids_leak'
    'ids_ipdc'
    'ids_command'
    'ids_limap'
    'ids_prio1'
    'ids_ldns'
    'ids_lirc'
    'ids_prio3'
    'ids_ipweb'
    'ids_pssh'
    'ids_lnetbios'
    'ids_misc'
    'ids_prio2'
    'ids_pnstd'
    'ids_lsql'
    'ids_psnmp'
    'ids_lssh'
    'ids_pnetbios3'
    'ids_lpop3'

Items in the folder:

- gpca.mat: data for the example
	- var_l: labels for the variables
	- obs_l: labels for the observations (timestamp when the interval ends)
	- weight: degree of relevance of the variables (published in journal)
	- weight_alt: degree of relevance of the variables (revised according to description in folder Data, more sensible)
	- x: data with 2350 observations on 95 variables

- run.m: Script to perform the EDA using the MEDA Toolbox commands. To run the script, the 
	current directory should be the one where the script is stored.

- fortreemap.mat: file to upload to http://nesg.ugr.es/meda-visualization/ to issue the treemap visualization

- Data: folder with information about original data (112 variables)

- html: web version of the example
 
