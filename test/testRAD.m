function tests = testRAD
    % Value tests for RAD.m against the shared reference values in expected.csv
    % Run with: runtests('test/testRAD.m')
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..'));
    lines = readlines(fullfile(here, 'series.csv'), 'EmptyLineRule', 'skip');
    series = struct();
    for i = 1:numel(lines)
        parts = split(lines(i), ',');
        series.(parts(1)) = str2double(parts(2:end));
    end
    testCase.TestData.series = series;
    testCase.TestData.expected = readtable(fullfile(here, 'expected.csv'), 'TextType', 'string');
end

function testMatchesReference(testCase)
    E = testCase.TestData.expected;
    warning('off', 'RAD:degenerate');
    cleanup = onCleanup(@() warning('on', 'RAD:degenerate'));
    for i = 1:height(E)
        f = RAD(testCase.TestData.series.(E.series(i)), E.tau(i), logical(E.doCentre(i)));
        if isnan(E.RAD(i))
            verifyTrue(testCase, isnan(f), sprintf('%s tau=%d', E.series(i), E.tau(i)));
        else
            verifyEqual(testCase, f, E.RAD(i), 'RelTol', 1e-12);
        end
    end
end

function testDefaults(testCase)
    x = testCase.TestData.series.ar1;
    verifyEqual(testCase, RAD(x), RAD(x, 1, true));
    verifyEqual(testCase, RAD(x), RAD(x'));
end

function testDegenerateWarns(testCase)
    verifyWarning(testCase, @() RAD(testCase.TestData.series.constant), 'RAD:degenerate');
end

function testOldArgumentOrderIsRejected(testCase)
    x = testCase.TestData.series.ar1;
    verifyError(testCase, @() RAD(x, 0), 'RAD:tau');
    verifyError(testCase, @() RAD(x, 1, 2), 'RAD:doCentre');
end
