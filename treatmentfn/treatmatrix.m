function m = treatmatrix(m, templatemat)
% Apply template and remove low intensity points
m = m.*templatemat;
m(m<0.1*max(m(:))) = 0; %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%

% Apply BWareaopen (matlab function)
m = m.*bwareaopen(m>0, freq2co(220)*time2co(0.1)); %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
end