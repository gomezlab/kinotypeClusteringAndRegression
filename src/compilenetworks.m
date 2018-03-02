function y = compilenetworks(finalnet, newnet)

%do not include self-interactions or repeat interactions (check node1 and
%node 2 since we are dealing with undirected interaction types)

for i = 1:length(newnet.node1)
    if sum( strcmp(newnet.node1(i),newnet.node2(i)) ) == 0
        %not a self-interaction

        if sum( strcmp(newnet.node1(i),unique([finalnet.node1,finalnet.node2])) ) == 0
            %new node1 is not in final node1 or final node2
            finalnet.node1(end+1) = newnet.node1(i);
            finalnet.node2(end+1) = newnet.node2(i);
        elseif sum( strcmp(newnet.node1(i),unique([finalnet.node1,finalnet.node2])) ) > 0
            %new node1 is in final node1 or node2
            if sum( strcmp(newnet.node2(i),unique([finalnet.node1,finalnet.node2])) ) == 0
                %new node2 is not in final node1 or final node2
                finalnet.node1(end+1) = newnet.node1(i);
                finalnet.node2(end+1) = newnet.node2(i);
            elseif sum( strcmp(newnet.node2(i),unique([finalnet.node1,finalnet.node2])) ) > 0
                %new node1 and new node2 are in final node1 or final node2
                %check each pairing
                %(newnode1,newnode2) to (finalnode1,finalnode2), done x1,y22
                %(newnode2,newnode1) to (finalnode1,finalnode2), done x2,y12
                %(newnode1,newnode2) to (finalnode2,finalnode1), done x3,y21
                %(newnode2,newnode1) to (finalnode2,finalnode1)  done x4,y11
                x1   = strcmp(newnet.node1(i),finalnet.node1);
                y22  = strcmp(newnet.node2(i),finalnet.node2(x1));
                
                x2   = strcmp(newnet.node2(i),finalnet.node1);
                y12  = strcmp(newnet.node1(i),finalnet.node2(x2));
                
                x3   = strcmp(newnet.node1(i),finalnet.node2);
                y21  = strcmp(newnet.node2(i),finalnet.node1(x3));
                
                x4   = strcmp(newnet.node2(i),finalnet.node2);
                y11  = strcmp(newnet.node1(i),finalnet.node1(x4));
                if sum(y22) == 0
                    if sum(y12) == 0
                        if sum(y21) == 0
                            if sum(y11) == 0
                                finalnet.node1(end+1) = newnet.node1(i);
                                finalnet.node2(end+1) = newnet.node2(i);
                            end
                        end
                    end
                end
                
            end
            
        end

    end
end

y = finalnet;

