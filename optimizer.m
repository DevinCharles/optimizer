function varargout = optimizer(varargin)
%
%     Copyright (C) 2015  Devin C Prescott
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%     
%     Author:
%     Devin C Prescott
%     devin.c.prescott@gmail.com
%% Usage
% [data,x,eqns]=optimizer('verbose',false,'update',true);
%% Data Shared Across Functions
    global data
    %% Input Data
    if nargin == 1
        fname = varargin{1};
        verbose = false;
    else
        p = inputParser();
        addOptional(p,'fname',strcat(pwd,'\','input_data.xlsx'));
        addOptional(p,'verbose',false);
        addOptional(p,'update',false);
        addOptional(p,'iterations',false);
        addOptional(p,'SOM',false);
        addOptional(p,'SOM_size',1e5);
        
        parse(p,varargin{:})
        
        fname = p.Results.fname;
        verbose = p.Results.verbose;
        update = p.Results.update;
        iter_on = p.Results.iterations;
        som = p.Results.SOM;
        som_sz = p.Results.SOM_size;
    end
    try
        rows = importData(fname);
    catch
        rows = importData();
    end
    
    %% Logical Arrays for Variables, Constants, Constraints, Equations
    var = cellfun(@(x) strcmpi(x,'var'),{data.type});
    con = cellfun(@(x) strcmpi(x,'con'),{data.type});
    cst = cellfun(@(x) strcmpi(x,'cst'),{data.type});
    aux = cellfun(@(x) strcmpi(x,'aux'),{data.type});
    obj = cellfun(@(x) strcmpi(x,'obj'),{data.type});
    eqn = cst|aux|obj;

    % Row Numbers of Variables to use to update the spreadsheet
    nums = 1:length({data.type});
    row_nums = nums(strcmpi({data.type},'var'))+4;
    %% Create Symbolic Variables and Constants
    syms(rows{logical(con+var)});
    %% Equations
    eqn_names = rows(eqn);
    % Create temporary symbols for equation names
    syms(rows{eqn});
    
    % Make a cell array of all equations
    eqns{length(eqn_names),1}=[];
    for i = 1:length(eqn_names)
        eqns{i} = eval(val(eqn_names{i}));
    end

    %% Create functions from symbolic equations
    % This replaces arbitrary symbolic variable names with x(1),x(2)...x(n) 
    % and replaces symbolic constants with their values in a function f(x)
    
    % Get indexes of Upper and Lower limits
    % We need to do this here before we convert strings to symbols
    ids_ll = ~strcmpi(cellfun(@num2str,{data.LL},...
        'UniformOutput',false),'NaN');
    ids_ul = ~strcmpi(cellfun(@num2str,{data.UL},...
        'UniformOutput',false),'NaN');
    % Get only limits for eqns
    %ids_ll = ids_ll(eqn);
    %ids_ul = ids_ul(eqn);
    
    var_names = {data(var).name};
%     con_names = {data(var).name};
    var_des = {data(var).description};
    
    % Substitute symbolic equation names with their actual equations
    for j = 1:length(eqns)
        eqns{j} = subs(eqns{j},eqn_names,eqns);
        eqns{j} = subs(eqns{j},{data(con).name},[data(con).value]);
        for i = 1:sum(var)
            % This is the variable we swap all other variables with: x(i)
            if ~som
                swap_var = strcat('x(',num2str(i),')');
            else
                swap_var = strcat('X',num2str(i));
            end
            eqns{j} = subs(eqns{j},var_names{i},swap_var);
            if verbose && j*i==1
                disp('-----------------')
                disp('Variable Swapping')
                disp('-----------------')
            end
            if verbose && j==1
                fprintf('Input Variable %s (%s) is now %s.\n',...
                    var_names{i},var_des{i},swap_var);
            end
        end
    end
    
        
    % Check in Lower and Upper Limits for Strings
    for j = 1:length({data.LL})
        if ischar(data(j).LL)
            % Change strins to their actual equations
            data(j).LL = subs(data(j).LL,eqn_names,eqns);
            % Replace constants with their values
            data(j).LL = subs(data(j).LL,{data(con).name},[data(con).value]);
            % Replace variables with the string x(j)
            for i = 1:sum(var)
                % This is the variable we swap all other variables with: x(i)
                if ~som
                    swap_var = strcat('x(',num2str(i),')');
                else
                    swap_var = strcat('X',num2str(i));
                end
                data(j).LL = char(subs(data(j).LL,var_names{i},swap_var));
            end
        end
        if ischar(data(j).UL)
            % Change strins to their actual equations
            data(j).UL = subs(data(j).UL,eqn_names,eqns);
            % Replace constants with their values
            data(j).UL = subs(data(j).UL,{data(con).name},[data(con).value]);
            % Replace variables with the string x(j)
            for i = 1:sum(var)
                % This is the variable we swap all other variables with: x(i)
                if ~som
                    swap_var = strcat('x(',num2str(i),')');
                else
                    swap_var = strcat('X',num2str(i));
                end
                data(j).UL = char(subs(data(j).UL,var_names{i},swap_var));
            end
        end
    end
    
    if verbose
        % Print x(i) substitutions in Equations
        disp('-----------------')
        disp('Equation Swapping')
        disp('-----------------')
        for j = 1:length(eqns)
            fprintf('%s = %s\n',eqn_names{j},char(eqns{j}))
        end
        
        % Print Changed UL & LL
        disp('-----------------')
        disp(' Limits Swapping ')
        disp('-----------------')
        for j = 1:length({data.LL})
            if ~isnumeric(data(j).LL)
                fprintf('%s = %s\n',strcat(data(j).name,' LL'),char(data(j).LL))
            end
            if ~isnumeric(data(j).UL)
                fprintf('%s = %s\n',strcat(data(j).name,' UL'),char(data(j).UL))
            end
        end
    end
    
    %% Self Organized Mapping 
    if som
        som_vars = arrayfun(@(ll,ul) linspace(ll,ul,floor((som_sz)^(1/sum(var)))),...
            [data(var).LL],[data(var).UL],'UniformOutput',false);
        
        % Replace *,/,^ with matrix style
        reps = {'*','/','^';'.*','./','.^'};
        fns = cellfun(@char,eqns,'UniformOutput',false);
        lls = cellfun(@num2str,{data.LL},'UniformOutput',false);
        uls = cellfun(@num2str,{data.UL},'UniformOutput',false);
        for i = 1:3
            fns = strrep(fns,reps{1,i},reps{2,i});
            lls = strrep(lls,reps{1,i},reps{2,i});
            uls = strrep(uls,reps{1,i},reps{2,i});            
        end
        eqns = fns;
        [data(:).LL] = lls{:};
        [data(:).UL] = uls{:};
        % Create ND grid
        [nd_vars{1:sum(var)}]=ndgrid(som_vars{:});
        %[X1,X2,X3,X4] = ndgrid(som_vars{:});
        % Unpack Variables
        for i = 1:sum(var)
            eval(strcat('X',num2str(i),'= nd_vars{',num2str(i),'};'));
        end
        % Evaluate all cst or obj functions at all values of x
        matrix = cellfun(@eval,eqns(cst(eqn)|obj(eqn)),'UniformOutput',false);
        lleqns = ids_ll(eqn); lleqns = lleqns(cst(eqn)|obj(eqn));
        uleqns = ids_ul(eqn); uleqns = uleqns(cst(eqn)|obj(eqn));
        % Remove values outside constraints
        ind1 = cellfun(@(a,b) a>=eval(num2str(b)),...
            matrix(lleqns),{data((obj|cst)&ids_ll).LL}',...
            'UniformOutput',false);
        ind2 = cellfun(@(a,b) a<=eval(num2str(b)),...
            matrix(uleqns),{data((obj|cst)&ids_ul).UL}',...
            'UniformOutput',false);
        ind = [ind1(:);ind2(:)];
        inds = true(size(ind{1}));

        for i = 1:length(ind)
%             disp(sum(reshape(ind{i},1,numel(ind{i}))));
            inds = ind{i}&inds;
%             disp(sum(reshape(inds,1,numel(inds))));
%             disp(' ');
        end
        if sum(reshape(inds,1,numel(inds)))<1
            disp('No options left.')
            disp('Are you attempting to maximize your objective by ')
            disp('(multiplying by -1) and it''s missing it''s bounds?')
            varargout{1}=[];
            return
        end
        % Add our variables to the matrix (we'll want to "see" these too)
        nvar = sum(var);
        neqn = sum(cst(eqn)|obj(eqn));
        for i = 1:nvar
            matrix{end+1}=nd_vars{i};
        end
        for i = 1:length(matrix)
            matrix{i}(~inds)= [];
        end
%         temp = cellfun(@(a) reshape(a,numel(nd_vars{1}),1),matrix,'UniformOutput',false);
        matrix = cell2mat(matrix)';
%         matrix = reshape(temp',length(temp)/(nvar+neqn),nvar+neqn);
        
        % Som Toolbox
        Icol = (1:numel(matrix(:,1)))';
        Struct.comp_names = [eqn_names(cst(eqn)|obj(eqn));var_names(:)];
        Struct.labels = cellstr(num2str(Icol));
        Struct.label_names = {'Matrix Index'};
        Struct.data = matrix;
        som_write_data(Struct,'SpringSOM.data');
        sD = som_read_data('SpringSOM.data');
        sD = som_normalize(sD,'var');
        mpsize = 'big';
        labeltype = 'add1';
        sM = som_make(sD,'mapsize',mpsize);
        sM = som_autolabel(sM,sD,labeltype);
        
        nmat = length(matrix);
        drawCompPlanes(sM,nvar,neqn,nmat)
        
        linkaxes(findobj('Type','Axes'),'xy');
        varargout{1}=matrix;
        return
    end
    
    %% Create function from string
    % Auxiliary, Constraint, & Objective Root Equations (no LL/UL applied)
    % These can be used later to display results
    % ie: Ax(1)+Bx(2) = 0 = f(x)
    % not Ax(1)+Bx(2) <= C 
    bse_funs = cellfun(@str2func,...
        strcat('@(x)',cellfun(@char,eqns,...
        'UniformOutput',false)),'UniformOutput',false);
        
    obj_funs = cellfun(@str2func,...
        strcat('@(x)',cellfun(@char,eqns(obj(eqn)),...
        'UniformOutput',false)),'UniformOutput',false);
    
    cst_funs_ll = cellfun(@str2func,...
        strcat('@(x)',...
        cellfun(@num2str,{data(cst&ids_ll).LL},'UniformOutput',false),...
        '-(',...
        cellfun(@char,eqns(cst(eqn)&ids_ll(eqn)),'UniformOutput',false)',...
        ')'),...
        'UniformOutput',false);
    
    cst_funs_ul = cellfun(@str2func,...
        strcat('@(x)',...
        cellfun(@char,eqns(cst(eqn)&ids_ul(eqn)),'UniformOutput',false),...
        '-(',...
        cellfun(@num2str,{data(cst&ids_ul).UL},'UniformOutput',false)',...
        ')'),...
        'UniformOutput',false);
    nlcsts = [cst_funs_ll(:);cst_funs_ul(:)];
    
    %% Initial Conditions
    x0 = val(var_names);
    %% Boundary Conditions
    lb = [data(var).LL];
    ub = [data(var).UL];
    %% Linear Constraints
    % TODO:
    % GET LINEAR FUNCTIONS HERE AS A MATRIX Ax<=b, RATHER THAN INSIDE THE
    % NONLINEAR FUNCTION AREA... MAYBE FASTER?
    %% NonLinear Constraints
    gfun = @(x) deal([...
        nlcsts{1}(x),...
        nlcsts{2}(x),...
        nlcsts{3}(x),...
        nlcsts{4}(x),...
        nlcsts{5}(x)...
        ],[]);
%     gfun = @(x) deal(...
%         nlcsts{:}(x),...
%         ,[]);
    %% Final Setup and Run
    if iter_on
        display = 'iter';
    else
        display = 'off';
    end
    options = optimoptions(...
        'fmincon',...
        'Algorithm','interior-point',...
        'Display',display,...
        'MaxFunEvals',500e3,...
        'MaxIter',100e3,...
        'TolCon',1e-3);


    % [x,fval,exitflag,output]=fmincon(f(x),x0,A,b,Aeq,beq,...
    %   lb,ub,nonlcon,options)
    [x,~,status,solver_output] = fmincon(obj_funs{1},x0,[],[],[],[],...
        lb,ub,gfun,options);

    %% Output
    if status ~= 1
        disp('!!!!!!!!!!!!!!!')
        disp('!Solver Failed!')
        disp('!!!!!!!!!!!!!!!')
        disp(solver_output.message)
    end
    if status >=0
        disp('')
        disp('=================')
        disp('Optimizer Results')
        disp('=================')
        disp('Design  Variables')
        disp('-----------------')
        for i = 1:length(x)
            fprintf('Design Variable %s (%s) is now %1.3e\n',...
                    var_names{i},var_des{i},x(i));
        end
        disp('-----------------')
        disp('Equation  Results')
        disp('-----------------')
        for i = 1:length(eqns)
            fprintf('%s = %1.3e\n',...
                    eqn_names{i},bse_funs{i}(x));
        end
        if update
           for i = 1:length(x)
               xlswrite(fname,x(i),'data',strcat('D',num2str(row_nums(i))));
           end
        end
    end
        
    varargout{1} = data;
    varargout{2} = x;
    varargout{3} = cell2struct(bse_funs,eqn_names,1);

end

function drawCompPlanes(sM,nvar,neqn,nmat)
    ntot = (neqn+nvar);
    % Variable Component Planes in 2x2 Layout
    for fig=1:ceil(ntot/4)
        figure
        if fig == ceil(ntot/4)
            fh(fig) = som_show(sM,'comp',4*fig-3:ntot,'norm','d');
            if nmat <= 50
                lab = som_show_add('label',sM,'subplot','all');
        %         set(lab,'Color','white')
            end
        else
            fh(fig) = som_show(sM,'comp',4*fig-3:4*fig,'norm','d');
            if nmat <= 50
                lab = som_show_add('label',sM,'subplot','all');
        %         set(lab,'Color','white')
            end
        end
    end
    nfig=fig;
    
%     % Equation Component Planes in 2x2 Layout
%     for fig=1:ceil(neqn/4)
%         figure
%         if fig == ceil(neqn/4)
%             fh(nfig+fig) = som_show(sM,'comp',nfig+4*fig-3:neqn,'norm','d');
%             if nmat <= 50
%                 lab = som_show_add('label',sM,'subplot','all');
%         %         set(lab,'Color','white')
%             end
%         else
%             fh(nfig+fig) = som_show(sM,'comp',nfig+4*fig-3:nfig+4*fig,'norm','d');
%             if nmat <= 50
%                 lab = som_show_add('label',sM,'subplot','all');
%         %         set(lab,'Color','white')
%             end
%         end
%     end
end


function rows = importData(varargin)
    global data
    if nargin == 1
        fname = varargin{1};
    else
    % Get the file
        [filename,pathname,~] = uigetfile('*.xlsx');
        fname = strcat(pathname,filename);
    end
    [~,~,args] = xlsread(fname,'data','C4:H100');
    del = cellfun(@(x) any(isnan(x(:))), args(:,1));
    args(del,:)=[];
    % Get the row and column names
    cols = args(1,:);
    rows = args(2:end,1);
    data = cell2struct(args(2:end,:),cols,2);
end

function varargout = swap(varargin)
    global data
    names = varargin{1};
    if ~iscell(names);names={names};end
    
    idx = zeros(1,length({data.name}));
    for i = 1:length(names)
        idx = logical(idx + cellfun(@(x) strcmpi(x,names{i}),{data.name}));
    end
    if nargin == 1
        varargout{1} = {data(idx).swap};
    elseif nargin == 2
        A = {varargin{2}};
        [data(idx).swap] = A{:}; 
    end
end

function varargout = val(varargin)
    global data
    names = varargin{1};
    if ~iscell(names);names={names};end
    
    idx = zeros(1,length({data.name}));
    for i = 1:length(names)
        idx = logical(idx + cellfun(@(x) strcmpi(x,names{i}),{data.name}));
    end
    if nargin == 1
        varargout{1} = [data(idx).value];
    elseif nargin == 2
        A = num2cell(varargin{2});
        [data(idx).value] = A{:}; 
    end
end

function varargout = LL(varargin)
    names = varargin{1};
    if ~iscell(names);names={names};end
    global data
    idx = zeros(1,length({data.name}));
    for i = 1:length(names)
        idx = logical(idx + cellfun(@(x) strcmpi(x,names{i}),{data.name}));
    end
    if nargin == 1
        varargout{1} = [data(idx).LL];
    elseif nargin == 2
        A = num2cell(varargin{2});
        [data(idx).LL] = A{:}; 
    end
end

function varargout = UL(varargin)
    names = varargin{1};
    if ~iscell(names);names={names};end
    global data
    idx = zeros(1,length({data.name}));
    for i = 1:length(names)
        idx = logical(idx + cellfun(@(x) strcmpi(x,names{i}),{data.name}));
    end
    if nargin == 1
        varargout{1} = [data(idx).UL];
    elseif nargin == 2
        A = num2cell(varargin{2});
        [data(idx).UL] = A{:}; 
    end
end