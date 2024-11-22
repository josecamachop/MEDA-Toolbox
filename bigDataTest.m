classdef bigDataTest < matlab.unittest.TestCase
    % Unit tests for the bigData functions

    methods(Test)
        function testiniLmodel(testCase)
            X = simuleMV(20,10,'LevelCorr',8);
            model = iniLmodel(X);
            ok = checkLmodel(model);
            assert((ok==true) )

            disp(model.lvs)
        end

        function testvarLpca(testCase)
            model = iniLmodel(simuleMV(20,10,'LevelCorr',8));
            model.lvs = 0:10;
            varLpca(model, 'Option', '0');
            varLpca(model, 'Option', '1');
            close all
        end

        function testvarLpls(testCase)
            X = simuleMV(20,10,'LevelCorr',8);
            Y = 0.1*randn(20,2) + X(:,1:2);
            model = iniLmodel(X,Y);
            model.lvs = 0:10;
            varLpls(model, 'Option', '0');
            varLpls(model, 'Option', '1');
            close all
        end

        function testupdateEwma(testCase)
        nobs = 100;
        nvars = 10;
        Lmodel = iniLmodel(simuleMV(nobs,nvars,'LevelCorr',6));
        Lmodel.type = 'PCA'; 
        Lmodel.prep = 2;  
        Lmodel.lvs = 1;
        Lmodel.nc = 100; % Number of clusters
        Lmodel.mat = loadingsLpca(Lmodel,'Option',0);
        mspcLpca(Lmodel);

        for i=1:4,
          nobst = 10;
          list(1).x = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(Lmodel.centr)*(nobst-1)/(Lmodel.N-1));
          Lmodel = updateEwma(list,'path',[],'Lmodel',Lmodel, 'debug', 0);
          mspcLpca(Lmodel);
        end
        close all
        end

        function testupdateIterative(testCase)
        nobs = 100;
        nvars = 10;
        Lmodel = iniLmodel;
        Lmodel.type = 'PCA'; 
        Lmodel.prep = 2;  
        Lmodel.lvs = 1;
        Lmodel.nc = 100; % Number of clusters
        Lmodel.path = '.\deleteMe\'; % MODIFY PATH TO YOUR CONVENIENCE

        for i=1:10,
          list(i).x = simuleMV(nobs,nvars,'LevelCorr',6);
        end

        Lmodel = updateIterative(list,'path',[],'Lmodel',Lmodel,'step',0.1,'files',1,'debug',0);
        mspcLpca(Lmodel);
        close all
        rmdir(Lmodel.path, 's')
        end

        function testleveragesLpca(testCase)
            X = simuleMV(20,10,'LevelCorr',8);
            Lmodel = iniLmodel(X);
            Lmodel.lvs = 1:3;
            L = leveragesLpca(Lmodel, 'Option', 0);
            L = leveragesLpca(Lmodel, 'Option', 1);
            close all
        end

        function testleveragesLpls(testCase)
        X = simuleMV(20,10,'LevelCorr',8);
        Y = 0.1*randn(20,2) + X(:,1:2);
        Lmodel = iniLmodel(X,Y);
        Lmodel.lvs = 1:3;
        L = leveragesLpls(Lmodel, 'Option', '0');
        L = leveragesLpls(Lmodel, 'Option', 1);
        close all
        end

        function testloadingsLpca(testCase)
        X = simuleMV(20,10,'LevelCorr',8);
        Lmodel = iniLmodel(X);
        Lmodel.lvs = 1:3;
        P = loadingsLpca(Lmodel, 'Option', 1, 'BlurIndex', .05);
        close all
        end

        function testloadingsLpls(testCase)
        X = simuleMV(20,10,'LevelCorr',8);
        Y = 0.1*randn(20,2) + X(:,1:2);
        Lmodel = iniLmodel(X,Y);
        Lmodel.lvs = 1;
        loadingsLpls(Lmodel);
        Lmodel.lvs = 1:2;
        loadingsLpls(Lmodel, 'Option', '0');
        close all
        end

        function testmedaLpca(testCase)
        X = simuleMV(20,10,'LevelCorr',8);
        Lmodel = iniLmodel(X);
        Lmodel.lvs = 1:3;
        map = medaLpca(Lmodel,'Threshold', 0.3,'Option','111');
        close all
        end

        function testmedaLpls(testCase)
        X = simuleMV(20,10,'LevelCorr',8);
        Y = 0.1*randn(20,2) + X(:,1:2);
        Lmodel = iniLmodel(X,Y);
        Lmodel.lvs = 1:3;
        map = medaLpls(Lmodel,'Threshold',0.3,'Option','111');
        close all
        end

        function testmspcAdicov(testCase)
        nobs = 100;
        nvars = 10;
        nPCs = 1;
        Lmodel = iniLmodel(simuleMV(nobs,nvars,'LevelCorr',6));
        Lmodel.multr = 100*rand(nobs,1); 
        Lmodel.lvs = 1:nPCs;

        nobst = 10;
        test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(Lmodel.centr)*(nobst-1)/(Lmodel.N-1));
        test(6:10,:) = 3*test(6:10,:);

        [Dstt,Qstt] = mspcAdicov(Lmodel,'Test', test(1:5,:),'Index',0);
        [Dstta,Qstta] = mspcAdicov(Lmodel,'Test',test(6:10,:),'Index',0);
        end

        function testmspcLpca(testCase)
        nobs = 100;
        nvars = 10;
        nPCs = 1;
        Lmodel = iniLmodel(simuleMV(nobs,nvars,'LevelCorr',6));
        Lmodel.multr = 100*rand(nobs,1); 
        Lmodel.lvs = 1:nPCs;

        nobst = 10;
        test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(Lmodel.centr)*(nobst-1)/(Lmodel.N-1));
        test(6:10,:) = 3*test(6:10,:);

        [Dst,Qst,Dstt,Qstt] = mspcLpca(Lmodel,'Test',test,'Option',100, ...
            'ObsLabel',[],'ObsClass',[1*ones(5,1);2*ones(5,1)]);
        close all
        end

        function testmspcLpls(testCase)
        nobs = 100;
        nvars = 10;
        nLVs = 1;
        X = simuleMV(nobs,nvars,'LevelCorr',6);
        Y = 0.1*randn(nobs,2) + X(:,1:2);
        Lmodel = iniLmodel(X,Y);
        Lmodel.multr = 100*rand(nobs,1); 
        Lmodel.lvs = 1:nLVs;

        nobst = 10;
        test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(Lmodel.centr)*(nobst-1)/(Lmodel.N-1));
        test(6:10,:) = 3*test(6:10,:);

        [Dst,Qst,Dstt,Qstt] = mspcLpls(Lmodel,'Test',test,'Option',100,'ObsClass',[1*ones(5,1);2*ones(5,1)]);

        close all
        end

        function testomedaLpca(testCase)
        nobs = 100;
        nvars = 10;
        nPCs = 10;
        Lmodel = iniLmodel(simuleMV(nobs,nvars,'LevelCorr',6));
        Lmodel.multr = 100*rand(nobs,1); 
        Lmodel.lvs = 1:nPCs;

        nobst = 10;
        test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(Lmodel.centr)*(nobst-1)/(Lmodel.N-1));
        test(1,1:2) = 10*max(abs(Lmodel.centr(:,1:2))); 
        dummy = zeros(10,1);
        dummy(1) = 1;

        omedavec = omedaLpca(Lmodel,test,dummy);
        close all
        end

        function testomedaLpls(testCase)
        nobs = 100;
        nvars = 10;
        nLVs = 10;
        X = simuleMV(nobs,nvars,'LevelCorr',6);
        Y = 0.1*randn(nobs,2) + X(:,1:2);
        Lmodel = iniLmodel(X,Y);
        Lmodel.multr = 100*rand(nobs,1); 
        Lmodel.lvs = 1:nLVs;

        nobst = 10;
        test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(Lmodel.centr)*(nobst-1)/(Lmodel.N-1));
        test(1,1:2) = 10*max(abs(Lmodel.centr(:,1:2))); 
        dummy = zeros(10,1);
        dummy(1) = 1;

        omedavec = omedaLpls(Lmodel,test,dummy);
        close all
        end

        function testscoresL(testCase)
        nobs = 100;
        nvars = 10;
        X = simuleMV(nobs,nvars,'LevelCorr',8);
        Lmodel = iniLmodel(X);

        nobst = 10;
        test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(X)*(nobst-1)/(nobs-1));

        Lmodel.multr = ceil(10*rand(nobs,1));

        Lmodel.lvs = 1;
        Lmodel = Lpca(Lmodel);
        scoresL(Lmodel,'Test',test);

        Lmodel.lvs = 1:2;
        Lmodel = Lpca(Lmodel);
        scoresL(Lmodel,'Test',test);
        close all
        end

        function testscoresLpca(testCase)
        nobs = 100;
        nvars = 10;
        X = simuleMV(nobs,nvars,'LevelCorr',8);
        Lmodel = iniLmodel(X);

        nobst = 10;
        test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(X)*(nobst-1)/(nobs-1));

        Lmodel.multr = ceil(10*rand(nobs,1));

        Lmodel.lvs = 1;
        scoresLpca(Lmodel,'Test',test);

        Lmodel.lvs = 1:2;
        [T,TT] = scoresLpca(Lmodel,'Test',test);
        % close all
        end

    end

end
