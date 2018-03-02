function [G] = ReadNetwork(nodisplay)

% by Scott Provan

% Requests and then reads a network file and puts it into the network structure form.
%
% INPUT FORMAT
%   NODELIST <list of node input fields>
%     <list of node labels together with appropriate node field parameters as specified>
%   one line of space
%   SOURCE-SINK <source and sink node names> (This line is optional.)
%   no line of space
%   EDGELIST <DIRECTED,UNDIRECTED> <list of arc input fields>
%     <list of appropriate arc field parameters as specified
%
% The entire network structure is displayed unless the nodisplay input is given any (nonnull) value

[fname,pname]=uigetfile('*.net','OPEN');
if fname==0, disp('No Network Read'), G=[]; return, end
fid=fopen(sprintf('%s%s',pname,fname),'r');
eval(sprintf('cd ''%s''',pname))
G=[];

line=1;

% Read node data

header=fgetl(fid);
parameter=sscanf(header,'%s',1);

if strcmp(parameter,'NETWORKTITLE')
    G.nettitle=header(13:end);
    header=fgetl(fid);
    parameter=sscanf(header,'%s',1);
end

if ~strcmp(parameter,'NODELIST')
    disp(sprintf('Improper input on Line %d: Looking for ''''NODELIST''''',line))
    fclose(fid); return
end

readflag=1; formatstring='%s'; formatadd='%*s ';
parameterorder=0;
while 1
    % Find out what node data is being input
    formatstring=sprintf('%s %s',formatadd, formatstring);
    parameter=sscanf(header,formatstring,1);
    if ~isempty(parameter)
        parameterorder=parameterorder+1;
        nodefields{parameterorder}=parameter;
    else
        break
    end
end

G.n=0;
readflag=1; namelength=0;
while readflag==1
    % Find out parameters for each node
    line=line+1;
    header=fgetl(fid);
    if isempty(sscanf(header,'%s'))
        if G.n==0
            disp(sprintf('Improper input on Line %d: No node input. (Did you put in a blank line?)',line))
            fclose(fid); return
        end
        readflag=0;
    else
        G.n=G.n+1;
        G.node(G.n).name=sscanf(header,'%s',1);
        namelength=max(namelength,length(G.node(G.n).name));
        formatstring='%*s'; formatadd='%*f';
        for i=1:parameterorder
            parameter=sscanf(header,[formatstring,' %f'],1);
            if isempty(parameter)
                disp(sprintf('Improper input on Line %d: Looking for parameter value',line))
                fclose(fid); return
            else
                eval(sprintf('G.node(G.n).%s = parameter;',nodefields{i}));
            end
            formatstring=sprintf('%s %s',formatstring,formatadd);
        end
    end
end

line=line+1;
header=fgetl(fid);

parameter=sscanf(header,'%s',1);

if strcmp(parameter,'SOURCE-SINK')
    
    source =sscanf(header,'%*s %s',1);
    sink =sscanf(header,'%*s %*s %s',1);
    if isempty(source) || isempty(sink)
        disp(sprintf('Improper input on Line %d: source or sink node not specified',line))
        fclose(fid); return
    end
    G.s=FindNode(G,source);
    G.t=FindNode(G,sink);
    if G.s==-1 || G.t==-1 
        disp(sprintf('Improper input on Line %d: source or sink node not found',line))
        fclose(fid); return
    end
   
    line=line+1;
    header=fgetl(fid);
    parameter=sscanf(header,'%s',1);
    
end

% Read arc data

if strcmp(parameter,'EDGELIST')

    parameter=sscanf(header,'%*s %s',1);
    switch parameter
        case 'DIRECTED'
            G.directed = 1;
        case 'UNDIRECTED'
            G.directed = 0;
        otherwise
            disp(sprintf('Improper input on Line %d: Looking for ''''DIRECTED'''' or ''''UNDIRECTED''''',line))
            fclose(fid); return
    end
    
    readflag=1; formatstring='%*s %s'; formatadd='%*s ';
    parameterorder=0;
 while 1
    % Find out what arc data is being input
    formatstring=sprintf('%s %s',formatadd, formatstring);
    parameter=sscanf(header,formatstring,1);
    if ~isempty(parameter)
        parameterorder=parameterorder+1;
        arcfields{parameterorder}=parameter;
    else
        break
    end
end
    
    G.m=0; G.arc=[];
    readflag=1;
    while readflag==1
        % Find out parameters for each arc
        line=line+1;
        header=fgetl(fid);
        if header==-1
            readflag=0;
        elseif isempty(sscanf(header,'%s'))
            readflag=0;
            if G.m==0
                disp(sprintf('Improper input on Line %d: No arc input. (Did you put in a blank line?)',line))
                fclose(fid); return
            end
        else
            G.m=G.m+1;
            %register head, tail
            tail=sscanf(header,'%s',1);
            head=sscanf(header,'%*s %s',1);
            if isempty(tail) | isempty(head)
                sprintf('Improper input on Line %d: Looking for arc endpoints',line)
                fclose(fid); return;
            end
            G.arc(end+1).tail=find(strcmp(tail,{G.node.name}));
            G.arc(end).head=find(strcmp(head,{G.node.name}));
            if length(G.arc(end).tail)~=1 | length(G.arc(end).tail)~=1
                sprintf('Improper input on Line %d: arc endpoints are either not on node list or duplicated in node list',line)
                fclose(fid); return;
            end
            
            formatstring='%*s %*s'; formatadd='%*f';
            for i=1:parameterorder
                parameter=sscanf(header,[formatstring,' %f'],1);
                if isempty(parameter)
                    sprintf('Improper input on Line %d: Looking for parameter value',line)
                    fclose(fid); return;
                else
                    eval(sprintf('G.arc(G.m).%s = parameter;',arcfields{i}));
                end
                formatstring=sprintf('%s %s',formatstring,formatadd);
            end
        end
       
    end
    
else
    
    disp(sprintf('Improper input on Line %d: Looking for ''''SOURCE-SINK'''' OR ''''EDGELIST''''',line))
    fclose(fid); return
    
end

fclose(fid); 
disp('Network successfully read')
if nargin==0, DisplayNetwork(G); end


