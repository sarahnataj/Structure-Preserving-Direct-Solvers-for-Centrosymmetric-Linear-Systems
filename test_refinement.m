clc
clear
ecase=2;
n=19:2:31;
for j=1:length(n)
[ e_d(j),e_s(j), e_r(j),n_ir(j), condition(j),condition_scaled(j),siz(j)]=Refinement_gmres_xy(n(j),ecase);
end
Results=table(n',siz', e_d', e_s',n_ir', e_r',condition',condition_scaled',...
        'VariableNames',{'number of nodes', 'size of the matrix', ...
        'error double','error single','#refinement iterations','error refinement',....
        'condition number','scaled condition number'})
semilogy(n,e_d,'k-*',n,e_s,'b-*',n, e_r, 'r-*','LineWidth',2)
if ecase==1
    title('1D biharmonic equation');
elseif ecase==2
    title('2D biharmonic equation');
elseif ecase==3
    title('2D variable-coefficient biharmonic equation');
end
xlabel('N')
ylabel('error')
legend('error double','error single','error refinement');
plotformat(1.5,6)