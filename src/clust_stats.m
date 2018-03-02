function [clust] = clust_stats(net, subnet, uniqdrugs, inhibunique90, n, underk)
% net          : formatted network
% d            : clusters (actual numbers)
% subnet       : name of subnetwork variable in "net"
% uniqdrugs    : drug list
% inhibunique90: drug targets 

d = net.subnets.(char(subnet));
%%% gather some cluster statistics
for i = 1:max(d)
    %determine how many kinases are in each cluster
    clust.clusters(i).name = strcat('cluster',int2str(i));
    clust.clusters(i).numkin = sum(d == i);
    %determine how many kinases are understudied in each cluster
    clust.clusters(i).numunderstud = length(intersect(n(d==i),underk));
    clust.clusters(i).understudcov = clust.clusters(i).numunderstud / length(n(d==i));
    
    clust.clusters(i).numdrugs = 0;
    
    %determine how many kinases are targeted in each cluster
    for j = 1:length(uniqdrugs)
        if ~isempty(intersect(n(d==i),inhibunique90.(char(uniqdrugs(j)))))
            clust.clusters(i).numdrugs = clust.clusters(i).numdrugs + 1;
        end
    end
    
end

%determine how many kinases in each cluster that do NOT have a drug
%targeting it
% for i = 1:length(clust.clusters)
%     x = n(d==i);
%     clust.clusters(i).numkintargeted = 0;
%     clust.clusters(i).numunderstudtarg = 0;
%     for j = 1:length(x)
%         y = FindNode(net,x(j));
%         if net.node(y).numdrugs > 0
%             clust.clusters(i).numkintargeted = clust.clusters(i).numkintargeted + 1;
%             if net.node(y).understudied == 1
%                 clust.clusters(i).numunderstudtarg = clust.clusters(i).numunderstudtarg + 1;
%             end
%         end
%     end
%     clust.clusters(i).drugcov = clust.clusters(i).numkintargeted / clust.clusters(i).numkin;
% end




% for k = 1:length(clust.clusters)
%     dchar = [];
%     dunder = [];
%     dtar = [];
%     duntar = [];
%     for i = 1:length(net.node)
%         if (net.node(i).(char(subnet)) == k) && (net.node(i).understudied == 0) 
%             dchar = [dchar; net.node(i).degree];
%         elseif (net.node(i).(char(subnet)) == k) && (net.node(i).understudied == 1)
%             dunder = [dunder; net.node(i).degree];
%         end
%         
%         if (net.node(i).(char(subnet))) == k && (net.node(i).numdrugs > 0)
%             dtar = [dtar; net.node(i).degree];
%         elseif (net.node(i).(char(subnet))) == k && (net.node(i).numdrugs == 0)
%             duntar = [duntar; net.node(i).degree];
%         end
%     end
%     
%     clust.clusters(k).medianchardegree = median(dchar);
%     clust.clusters(k).medianunderdegree = median(dunder);
%     clust.clusters(k).mediantargetdegree = median(dtar);
%     clust.clusters(k).medianuntargetdegree = median(duntar);
%     
% end
