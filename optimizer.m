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
        parse(p,varargin{:})
        fname = p.Results.fname;
        verbose = p.Results.verbose;
        update = p.Results.update;
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
    eqn = logical(cst+aux+obj);
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
            swap_var = strcat('x(',num2str(i),')');
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
                swap_var = strcat('x(',num2str(i),')');
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
                swap_var = strcat('x(',num2str(i),')');
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
        cellfun(@char,eqns(cst(eqn)&ids_ll(eqn)),'UniformOutput',false),...
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
    if verbose
        display = 'iter';
    else
        display = 'off';
    end
    options = optimoptions(...
        'fmincon','Algorithm',...
        'interior-point',...
        'Display',display,...
        'MaxFunEvals',100000,...
        'MaxIter',30000,...
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
            fprintf('Design Variable %s (%s) is now %s\n',...
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