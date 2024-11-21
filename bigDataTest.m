classdef bigDataTest < matlab.unittest.TestCase
    % Unit tests for the bigData functions

    methods(Test)
        function testiniLmodel(testCase)
            X = rand(20, 10);
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
        Lmodel.mat = loadingsLpca(Lmodel,0);
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
    end

end