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
            X = rand(20, 10);
            model = iniLmodel(X);
            model.lvs = 1:2;
            varLpca(model, 'Option', '0')
            varLpca(model, 'Option', '1')
        end
    end

end
