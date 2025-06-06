classdef plotVecTest < matlab.unittest.TestCase
    % Unit tests for the plotVec test function

    methods(Test)
        function testplotVec(testCase)
            X = simuleMV(26,1,'LevelCorr',8);
            label = 'A':'Z';
            label = label(1:26); 
            
            % Categorical Classes
            class_cat = repmat({'consonant'}, 26, 1);
            vocal_id = [1, 5, 9, 15, 21];
            class_cat(vocal_id) = {'vocal'};
            class_cat{25} = '?'; 
            
            % Numerical (continuous) Classes
            class_num = ((1:26)/2)';
            class_num(1) = 20; % Letter 'A' should dominate the colors

            % Multiplicity
            mult = [repmat(.5,1,4), ...
                    repmat(15,1,4), ...
                    repmat(30,1,4), ...
                    repmat(75,1,4), ...
                    repmat(125,1,4), ...
                    repmat(200,1,10), ...
                    200 ];

            % Alpha
            alphas = [repmat(.3,1,4), ...
                        repmat(.5,1,4), ...
                        repmat(.7,1,4), ...
                        repmat(.9,1,4), ...
                        repmat(1,1,4), ...
                        repmat(.1,1,10)];

            % Categorical classes
                plotVec(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_cat);
                title("Categorical classes - Okabe as default color if the number of classes is <=8", "FontSize", 10)

            % Line plot - Categorical
                plotVec(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_cat, 'PlotType', 'Lines');
                title("Categorical classes - Okabe as default color if the number of classes is <=8", "FontSize", 10)
            
            % Categorical classes - Alpha
                plotVec(X, 'XYLabel', {"X", "Y"}, ...
                'EleLabel', label, 'ObsClass', class_cat, 'ObsAlpha', alphas(1:26));
                title("Categorical classes - Opacity", "FontSize", 10)

            % Numerical (continuous) classes
                plotVec(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num);
                title("Numerical continuous classes - Non equidistant colors. A = 20", "FontSize", 10)

            % Line plot - Numerical
            plotVec(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num, 'PlotType', 'Lines');
                title("Numerical continuous classes - Non equidistant colors. A = 20", "FontSize", 10)

            % Multiplicity  
                plotVec(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num, 'Multiplicity', mult(1:26));
                title("Numerical continuous classes - Multiplicity", "FontSize", 10)

            % Numerical classes - Alpha
                plotVec(X, 'XYLabel', {"X", "Y"}, 'EleLabel', label, ...
                'ObsClass', class_num, 'Multiplicity', mult(1:26), 'ObsAlpha', alphas(1:26));
                title("Numerical continuous classes - Opacity", "FontSize", 10)
        end


        function testplotVecM(testcase)
            % plotVec tests for input of shape (N, M), with M > 1

            % bars/lines with constant control limits
            plotVec(randn(100,3),'XYLabel',{'Functions','Time'},'LimCont',[1, -1, 3],'Color','parula');
            plotVec(randn(100,3),'XYLabel',{'Functions','Time'},'LimCont',[1, -1, 3],'Color','parula','PlotType','Lines');
            
            % many bars/lines with constant control limits
            plotVec(randn(100,30),'XYLabel',{'Functions','Time'},'LimCont',[1, -1, 3],'Color','parula');
            plotVec(randn(100,30),'XYLabel',{'Functions','Time'},'LimCont',[1, -1, 3],'Color','parula','PlotType','Lines');

            % bars with different opacity values
            N = 10;
            alphas = rand(N, 1);
            alphas(fix(N/3):fix(2*N/3)) = 1;
            plotVec(randn(N,3),'XYLabel',{'Functions','Time'},'LimCont',1,'Color','summer','ObsAlpha',alphas);
            
            % with labels and categorical classes in observations and variable limit
            plotVec(randn(5,3),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{[],'Functions'},'LimCont',randn(5,1));
            plotVec(randn(5,3),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{[],'Functions'},'LimCont',randn(5,1),'PlotType','Lines');           
            
            % numerical classes in observations and variable limit           
            plotVec(randn(10,3),'ObsClass',1:10,'ClassType','Numerical','XYLabel',{[],'Functions'},'LimCont',randn(10,1));
            plotVec(randn(10,3),'ObsClass',1:10,'ClassType','Numerical','XYLabel',{[],'Functions'},'LimCont',randn(10,1),'PlotType','Lines');
            
            % labels, multiplicity and classes in observations and variable limit
            plotVec(randn(5,3),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{[],'Functions'},'LimCont',randn(5,1),'Multiplicity',100*rand(5,1),'Markers',[20 50 100]);
            plotVec(randn(5,3),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{[],'Functions'},'LimCont',randn(5,1),'Multiplicity',100*rand(5,1),'Markers',[20 50 100],'PlotType','Lines');
        end

        
        function testvarPca(testCase)
            X = simuleMV(26,30,'LevelCorr',8);
            label = 'A':'Z';
            label = label(1:26); 
            varPca(X)
        end
        
    end
end