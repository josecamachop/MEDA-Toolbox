classdef dataVizTest < matlab.unittest.TestCase
    % Unit tests for the bigData functions

    methods(Test)
        function testplotScatter(testCase)
            X = simuleMV(26,2,'LevelCorr',8);
            label = 'A':'Z';
            label = label(1:26); 
            
            % Categorical Classes
            class_cat = repmat("consonant", 26, 1);
            vocal_id = [1, 5, 9, 15, 21];
            class_cat(vocal_id, :) = "vocal";
            class_cat(25, :) = "?";

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

            % Categorical classes
                % Okabe default color (classes <=8)  
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_cat);
                title("Categorial classes - Okabe as default color if the number of classes is <=8", "FontSize", 10)

            % Numerical classes
                % parula default color   
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num);
                title("Numerical continuous classes - Non equdistant colors. A = 20", "FontSize", 10)
                
            % Plot type: "size"
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num,...
                'PlotMult', 'size', 'Multiplicity', mult(1:26));
                title("Numerical continuous classes - Multiplicity as size", "FontSize", 10)
            
            % Plot type: "shape"
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num,...
                'PlotMult', 'shape', 'Multiplicity', mult(1:26));
                title("Numerical continuous classes - Multiplicity as marker shape", "FontSize", 10)
            
            
            % Plot type: "zaxis"
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num,...
                'PlotMult', 'zaxis', 'Multiplicity', mult(1:26));
                title("Numerical continuous classes - Multiplicity in Z axis", "FontSize", 10)
            
            % Plot type: "zsize"
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num,...
                'PlotMult', 'zsize', 'Multiplicity', mult(1:26));
                title("Numerical continuous classes - Multiplicity in Z axis and shape.", "FontSize", 10)
            
            % Plot type: "zshape"
                plotScatter(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num,...
                'PlotMult', 'zshape', 'Multiplicity', mult(1:26));
                title("Numerical continuous classes - Multiplicity in Z axis and shape.", "FontSize", 10)
        end

        function testplotVec(testCase)
            X = simuleMV(26,1,'LevelCorr',8);
            label = 'A':'Z';
            label = label(1:26); 
            
            % Categorical Classes
            class_cat = repmat("consonant", 26, 1);
            vocal_id = [1, 5, 9, 15, 21];
            class_cat(vocal_id, :) = "vocal";
            class_cat(25, :) = "?";
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

            % Categorical classes
                plotVec(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_cat);
                title("Categorial classes - Okabe as default color if the number of classes is <=8", "FontSize", 10)
 
            % Line plot - Categorical
                plotVec(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_cat, 'PlotType', 'Lines');
                title("Categorial classes - Okabe as default color if the number of classes is <=8", "FontSize", 10)

            % Numerical (continuous) classes
                plotVec(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num);
                title("Numerical continuous classes - Non equdistant colors. A = 20", "FontSize", 10)

            % Line plot - Numerical
            plotVec(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num, 'PlotType', 'Lines');
                title("Numerical continuous classes - Non equdistant colors. A = 20", "FontSize", 10)

            % Multiplicity  
                plotVec(X, 'XYLabel', {"X", "Y"}, ... 
                'EleLabel', label, 'ObsClass', class_num, 'Multiplicity', mult(1:26));
                title("Numerical continuous classes - Multiplicity", "FontSize", 10)

        end


        function testvarPca(testCase)
            X = simuleMV(26,30,'LevelCorr',8);
            label = 'A':'Z';
            label = label(1:26); 
            varPca(X)
        end

    end

end

