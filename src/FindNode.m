function [index] = FindNode(G,nodename)

% Finds index of node with name nodename, -1 if no such node

index = find(strcmp(nodename,{G.node.name}));
if isempty(index), index=-1; end