clc
clear
ecase=2;
n=19:2:30;
for j=1:length(n)
[ e_ub(j), e_u(j),n_ir(j), condition(j),condition_scaled(j),siz(j)]=Refinement_gmres_xy(n(j),ecase);
end
Results=table(n',siz', e_ub',n_ir', e_u',condition',condition_scaled',...
        'VariableNames',{'number of nodes', 'size of the matrix', ...
        'error A\b','#refinement iterations','error refinement',....
        'condition number','scaled condition number'})
semilogy(n,e_ub,'-*',n, e_u, '-*','LineWidth',2)
if ecase==1
    title('1D biharmonic equation');
elseif ecase==2
    title('2D biharmonic equation');
elseif ecase==3
    title('2D variable-coefficient biharmonic equation');
end
xlabel('N')
ylabel('error')
legend('error Ab','error refinement');
plotformat(1.5,6)