function do(thecmd)

disp(thecmd)
try
    eval(thecmd)
catch
    disp('FAILURE: ')
    disp(thecmd);
end