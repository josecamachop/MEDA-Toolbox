load SwA.ifInOctets.mod
load SwB.ifInOctets.mod
load SwC.ifInOctets.mod
load SwA.ifOutOctets.mod
load SwB.ifOutOctets.mod
load SwC.ifOutOctets.mod

Todo = [SwA_ifInOctets(end-103:end,:) SwA_ifOutOctets(end-103:end,:) SwB_ifInOctets(end-103:end,:) SwB_ifOutOctets(end-103:end,:) SwC_ifInOctets(end-103:end,:) SwC_ifOutOctets(end-103:end,:)];
Todo = Todo(2:104,:)-Todo(1:103,:);
lab={'SwA.ifIn1', 'SwA.ifIn2', 'SwA.ifIn3', 'SwA.ifIn4', 'SwA.ifIn5', 'SwA.ifIn6', 'SwA.ifIn7', 'SwA.ifIn8', 'SwA.ifIn9', 'SwA.ifIn10', 'SwA.ifIn11', 'SwA.ifIn12', 'SwA.ifIn13', 'SwA.ifIn14', 'SwA.ifOut1', 'SwA.ifOut2', 'SwA.ifOut3', 'SwA.ifOut4', 'SwA.ifOut5', 'SwA.ifOut6', 'SwA.ifOut7', 'SwA.ifOut8', 'SwA.ifOut9', 'SwA.ifOut10', 'SwA.ifOut11', 'SwA.ifOut12', 'SwA.ifOut13', 'SwA.ifOut14', 'SwB.ifIn1', 'SwB.ifIn2', 'SwB.ifIn3', 'SwB.ifIn4', 'SwB.ifIn5', 'SwB.ifIn6', 'SwB.ifIn7', 'SwB.ifIn8', 'SwB.ifIn9', 'SwB.ifIn10', 'SwB.ifIn11', 'SwB.ifIn12', 'SwB.ifIn13', 'SwB.ifIn14', 'SwB.ifOut1', 'SwB.ifOut2', 'SwB.ifOut3', 'SwB.ifOut4', 'SwB.ifOut5', 'SwB.ifOut6', 'SwB.ifOut7', 'SwB.ifOut8', 'SwB.ifOut9', 'SwB.ifOut10', 'SwB.ifOut11', 'SwB.ifOut12', 'SwB.ifOut13', 'SwB.ifOut14', 'SwC.ifIn1', 'SwC.ifIn2', 'SwC.ifIn3', 'SwC.ifIn4', 'SwC.ifIn5', 'SwC.ifIn6', 'SwC.ifIn7', 'SwC.ifIn8', 'SwC.ifIn9', 'SwC.ifIn10', 'SwC.ifIn11', 'SwC.ifIn12', 'SwC.ifIn13', 'SwC.ifIn14', 'SwC.ifOut1', 'SwC.ifOut2', 'SwC.ifOut3', 'SwC.ifOut4', 'SwC.ifOut5', 'SwC.ifOut6', 'SwC.ifOut7', 'SwC.ifOut8', 'SwC.ifOut9', 'SwC.ifOut10', 'SwC.ifOut11', 'SwC.ifOut12', 'SwC.ifOut13', 'SwC.ifOut14'};

Todo = Todo(2:end,:);
cal = Todo([1:3 5:38 92:end],:);
test = Todo([39:91],:);

i_d = find(sum(cal,1)>0);
i_d = i_d([1:33 35:end]);
cal = cal(:,i_d);
test = test(:,i_d);
lab = lab(i_d);

save datos_proc cal test lab