classdef plotScatterTest < matlab.unittest.TestCase
    % Unit tests for the plotScatter function

    methods(Test)
        function testplotScatter(testCase)
            X = simuleMV(26,2,'LevelCorr',8);
            label = 'A':'Z';
            label = label(1:26); 
            
            % Categorical Classes
            class_cat = repmat({'consonant'}, 26, 1);
            vocal_id = [1, 5, 9, 15, 21];
            class_cat(vocal_id) = {'vocal'};
            class_cat{25} = '?'; 

            % Numerical (continuous) Classes
            class_num = (((1:26)/2))';
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
                % Okabe default color (classes <=8)  
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_cat, 'ObsAlpha', alphas(1:26));
                title("Categorical classes - Okabe as default color if the number of classes is <=8", "FontSize", 10)

            % Numerical classes
                % parula default color   
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num, 'ObsAlpha', alphas(1:26));
                title("Numerical continuous classes - Non equidistant colors. A = 20", "FontSize", 10)
                
            % Plot type: "size"
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num,...
                'PlotMult', 'size', 'Multiplicity', mult(1:26), 'ObsAlpha', alphas(1:26));
                title("Numerical continuous classes - Multiplicity as size", "FontSize", 10)
            
            % Plot type: "shape"
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num,...
                'PlotMult', 'shape', 'Multiplicity', mult(1:26), 'ObsAlpha', alphas(1:26));
                title("Numerical continuous classes - Multiplicity as marker shape", "FontSize", 10)
            
            
            % Plot type: "zaxis"
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num,...
                'PlotMult', 'zaxis', 'Multiplicity', mult(1:26), 'ObsAlpha', alphas(1:26));
                title("Numerical continuous classes - Multiplicity in Z axis", "FontSize", 10)
            
            % Plot type: "zsize"
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num,...
                'PlotMult', 'zsize', 'Multiplicity', mult(1:26), 'ObsAlpha', alphas(1:26));
                title("Numerical continuous classes - Multiplicity in Z axis and shape.", "FontSize", 10)
            
            % Plot type: "zshape"
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num,...
                'PlotMult', 'zshape', 'Multiplicity', mult(1:26), 'ObsAlpha', alphas(1:26));
                title("Numerical continuous classes - Multiplicity in Z axis and shape.", "FontSize", 10)
        end

        function testplotScatter2(testCase)
            % random data with filled marks and control limits
            plotScatter(rand(100,2),'XYLabel',{'Y','X'},'LimCont',{0.8,0.8});
            
            % labels and classes in elements
            plotScatter(randn(5,2),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'},'Color','hsv');
            
            % labels, and numerical and categorical classes
            X = randn(5,2);
            plotScatter(X,'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'},'ClassType','Categorical');
            plotScatter(X,'EleLabel',{'one','two','three','four','five'},'ObsClass',1:5,'XYLabel',{'Y','X'},'ClassType','Numerical');
            
            % labels, multiplicity and classes in elements
            X = randn(5,2);
            plotScatter(X,'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'},'ClassType','Categorical','Multiplicity',[1 20 50 100 1000]);
            plotScatter(X,'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'},'ClassType','Numerical','Multiplicity',[1 20 50 100 1000]);
            
            % different plot types for multiplicity
            mult = {'size','shape','zaxis','zsize'};
            for o = 1:length(mult)
                plotScatter(X,'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'},'PlotMult',mult{o},'Multiplicity',[1 20 50 100 1000]);
            end
        end

    end
end

