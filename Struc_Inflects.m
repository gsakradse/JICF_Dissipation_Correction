function [Sdiss, plotStr] = Struc_Inflects(cS, dx)



cPiv = cS;
Mdx = dx;

% Find inflection points of the compensated structure functions
kk = 0;
for ii = 1:length(cS(1,:))
    
    logi = dx(:,ii) >= 0.1005 & dx(:,ii) <= 7;    
    
    offSet = find(logi > 0, 1);
    
    dCs(:,ii) = diff(cS(:,ii));
    
    % If the first derivative does not cross 0 in range of interest(
    % <=10-1m)
    if isempty(find(dCs(logi,ii) <= 0, 1))
        
        kk = kk+1;
        % Take the second derivative 
        
        ddCs = diff(dCs(logi,ii));
        int(ii) = find(ddCs >=0, 1) + offSet;
        
        fprintf('\nyo yo yo watch yourself!\nderivative of the strucute function in index %d does not cross zero', ii)
        fprintf('\nusing the first zero crossing of the 2nd derivative instead\n')
    else
        
    int(ii) = find(dCs(logi,ii) <= 0, 1) + offSet;
    
    end
    
    Sdiss(ii) = cS(int(ii),ii);
    
end


plotStr.Ndx = dx;
plotStr.cS = cS;
plotStr.dCs = dCs;
plotStr.int = int;


%% plot some things 
% 


%{
inds = [1 5 10 15 16 17 18 20 25 30];

cols = pc3(length(inds));

fpos = FigPosition([700 500], 1);
figure('pos', fpos)

kk = 1;
for ii = inds
   
    semilogx(Mdx(1:length(dCs),ii), dCs(:,ii), '-','color', cols(kk,:), ...
        'DisplayName', num2str(ii),'LineWidth', 2, 'MarkerSize', 5 ) 
    hold on
    semilogx(Mdx(int(ii), ii), dCs(int(ii), ii), 'ro', 'MarkerSize', 5)
   kk = kk+1; 
end
hold off
title('d/dr(Compensated Structure Function)', 'interpreter', 'latex')
legend show
lgnd = findobj(gcf, 'Tag', 'legend');
lgnd.Location = 'best';
lgnd.NumColumns = 3;


figure('pos', fpos)
kk = 1;
for ii = inds
   
    semilogx(Mdx(:,ii), cS(:,ii), '-','color', cols(kk,:),...
        'DisplayName', num2str(ii), 'LineWidth', 2, 'MarkerSize', 5 )
    hold on
    semilogx(Mdx(int(ii), ii), cS(int(ii), ii), 'ro', 'MarkerSize', 5)
    kk = kk+1;
end
hold off
title('Compensated Structure Function', 'interpreter', 'latex')
legend show
lgnd = findobj(gcf, 'Tag', 'legend');
lgnd.Location = 'best';
lgnd.NumColumns = 3;
%}