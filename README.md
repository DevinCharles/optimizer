# optimizer

Optimizer is a wrapper function for *fmincon* and the *[SOM-Toolbox](https://github.com/ilarinieminen/SOM-Toolbox).* It's primary goal is to make optimization easily accesible and intuitive for everyone. 

The function gets its all its input from a simple excel spreadsheet where you layout your design variables with initial values and bounds, constants, auxilary equations (equations that might make up your constraints or objective function), your constraint functions, and your objective function. There is a blank template included, and a sample file to show how things work.

If you're human, and are unsure where to start your initial values, there is a handy *discrete* mode, which will solve all your constraint equations and your objective function over your bounds over a discrete matrix, so you can get an idea of where to start.

If you're still human, and staring at a large matrix of numbers leaves you stumped, you can make use of the Self-Organizing Maps toolbox written by Esa Alhoniemi, Johan Himberg, Jukka Parviainen and Juha Vesanto. This will visually direct you to a good set of initial conditions, or even a "good enough" solution.

### Quick Start - Excel Input Sheet Example
| name|    value   | type | LL | UL | description |
|:---:|:----------:|:----:|:--:|:--:|:------------|
| A   |4           | con  |    |    | Some Constant Value |
| B   |3           | var  | 1  | 6  | A Design Variable Value with lower and upper bounds |
| C   |2           | var  |    |10  | A Design Variable with only one bound|
| D2  | C+B        | cst  |    |2*A | A Constraint Equation with an equation as a bound|
| Ei  | B-C^2      | cst  | -1 |    | An additional Constraint Equation|
| Fo  | A^2-D      | aux  |    |    | An auxilary equation (these should not have bounds)|
| Gm  | A*Fo-B^2+D2| obj  |  1 |    | The objective function**
** NOTE: The bound(s) on the objective function is(are) optional in *discrete* mode, this can help you narrow in on your solution. It will be ignored in *optimize* mode. Also, objective functions are always **minimized** in *optimize* mode, so if you need it maximized, multiply the whole equation by -1. Be careful to remove that for *discrete mode* if you have a bound, or make your bounds negative as well. I'm working on a fix... 

### Quick Start - Matlab Code
This is the quickest way to get started. It will run in *discrete* mode and give you a table of possible values for your solution. 

`[Table] = optimizer();` 

To run the *fmincon* optimization on your problem, use the below command. The *data* output is your input data from the excel spreadsheet in a Structure, x is the actual optimized design variables, and eqn is a structure with all your constraint, auxilary, and objective functions, as matlab functions of the design variables x(1), x(2)...

`[data,x,eqn] = optimizer();`

### Name-Value Input Options
- **Name [Default Value] - (Alternate Values)**
- *fname* [`'input_data.xlsx'`] - (`'C:\Some Full Path\file.xlsx'`)
  - This supplies the full path to your input excel file.
- *verbose* [`false`] - (`true`) 
  - This turns on and off verbose mode, which will show you your equation substitutions and will show you which design varible bounds are removing the most solutions.
- *update* [`false`] - (`true`) 
  - In *discrete* mode, *update* will create the matrix of unique solutions found (values to constraint and objective functions) to sheet 2 in your input data sheet. 
  - In *optimize* mode, *update* will update the initial conditions in your input data sheet with the optimal values it finds.
- *iterations* [`false`] - (`true`) 
  - This will show each iteration in *optimize* mode 
- *mode* [`'discrete'`] - (`'optimize'`) 
  - This is where you choose whether your want to use the discrete or optimize modes. Discrete is helpful for finding initial conditions "by hand" while *optimize* will "do it for you" to the best of it's ability.
- *disc_size* [`1e5`] - (`Any Integer`) 
  - This is the maximum size of your discrete matrix. r<sup>n</sup> = disc_size where *r* is the range between your bounds for each variable, and n is the number of variables. The greater this number, the more discrete values between your variable bounds, but the more memory intensive the operation. I suggest not exceeding 1e7. 
- *SOM* [`false`] - (`true`)
  - This will use the SOM-Toolbox to show component planes with your discrete solutions self organized. Basically, that means it will make color maps for each constraint, your objective, and each design variable. Each unique solution is mapped to the same position on each plane or map, and the color shows you its value in that plane. I've found this incredibly useful for finding solid initial conditions or for "hand picking" optimal values.

### Example
`[Table] = optimizer('fname',fname,'SOM',true,'disc_size',1e6,'verbose',true,'update',true);`
