%% Author: Li Muchen / Andrew Elizabeth
%% Date: 2026-03-18
%% Version: 1.0
%% License: MIT License

% Copyright (c) 2026 Li Muchen / Andrew Elizabeth / EasterGhost

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%% Before start
% For those who download this code, I warmly welcome you to study and use it!
% However, please be sure to comply with the MIT License, retain the copyright notice and license information, and refrain from deleting or modifying these details.
% I hope you can download, use, and disseminate this code in a standardized manner, respecting the original author's labor achievements.
% At the same time, we welcome your improvements and optimizations to the code, and you are encouraged to share your improved versions under the premise of complying with the license.
% If you have any questions, suggestions, or improvement ideas regarding the code, feel free to submit an issue. I will do my best to respond and assist you.
% Additionally, you are welcome to contribute your improved version by submitting a pull request on GitHub.
% Finally, thank you for your attention and support to this code! If you like this code, feel free to give it a star!
% I hope it can help you better understand and apply the relevant knowledge and techniques of structural topology optimization.

%% ---- Topology Optimization using Optimality Criteria Method ----
% This code is a modified version of the original top99.m, which is a MATLAB implementation of the topology optimization problem.
% This optimized code uses vectorized operations and persistent variables to optimize the topology of a structure.
% By using indexed access and pre-computed matrices(e.g. KE->K), it improves performance.
% This code's core idea of optimization is similar to the top88.m, which is vectorization and persistent variables to improve multi-core performance and avoid unnecessary recalculations.
% But it is easier to understand and more efficient than the top88.m in small and medium scale problems(<100k elements, large scale not tested yet).
% BTW, I really don't know why the original top88.m is a little bit slower.

%% ---- Introduction of Parameters ----
% nelx: Number of elements in the x-direction
% nely: Number of elements in the y-direction
% volfrac: Volume fraction of the material
% penal: Penalization factor for the design variables
% rmin: Filter radius for the sensitivity analysis
% nelem: Total number of elements in the mesh, computed as nelx * nely
% x: Design variable matrix, initialized with the volume fraction, which represents the pseudo-density of each element
% x_old: Previous design variable matrix, used to calculate the change in design variables
% x_new: Updated design variable matrix, computed using the Optimality Criteria method
% F: Force vector, initialized with a unit load at the center of the right edge to simulate a simple loading condition (Cantilever Beam)
% KE: Element stiffness matrix, defined for a 2D linear elastic element
% U: Displacement vector, initialized to zero
% K: Global stiffness matrix, initialized as a sparse matrix
% fixeddofs: Fixed degrees of freedom, defined for a simple fixed boundary condition where the left edge is fixed and the right edge is loaded
% freedofs: Free degrees of freedom, defined as the set difference between all degrees of freedom
% alldofs: All degrees of freedom, defined for the mesh
% edof: Element degrees of freedom, defined for each element in the mesh
% c: Objective function value, initialized to zero
% dc: Sensitivity of the objective function with respect to the design variables, initialized to zero
% dc_vec: Vectorized sensitivity of the objective function with respect to the design variables, computed as -penal * (x_vec .^(penal - 1)) .* UKU_all
% loop: Iteration counter, initialized to zero
% change: Change in design variables, initialized to 1.0 to enter the optimization loop
% H: Filter matrix, used to average the sensitivities over a neighborhood defined by rmin
% Hs: Row sums of the filter matrix, used for numerical stability
% iK, jK: Indices for the global stiffness matrix, pre-computed for efficiency
% iH, jH, sH: Indices and values for the filter matrix, pre-computed for efficiency
% idx(used in FE to build indices iK, jK for K): Index for the element degrees of freedom, begin at 1
% idx(used in sensitivity filtering to build indices iH, jH for H): Index for the sensitivity calculation, begin at 1
% Ue_all: Displacements for all elements, extracted by edof
% UKU_all: Element strain energy, computed as U'KU for all elements
% dcn: Filtered sensitivity, computed using the filter matrix H
% dcn_vec: Vectorized filtered sensitivity, computed as (H * (x_vec .* dc_vec)) ./ (Hs .* x_vec)
% x_vec: Vectorized design variable matrix, computed as x(:)
% x_penal: Penalized design variable matrix, computed as x .^ penal
% x_vec: Vectorized design variable matrix, computed as x(:)
% KE_vec: Vectorized element stiffness matrix, computed as KE(:)
% sK: Sparse global stiffness matrix, computed using KE(:) and x_penal
% nfilter: Total number of filter entries, computed as nelem * (2*floor(rmin) + 1)^2

%% ---- Main function to run the topology optimization problem ----
function top99_optimized(nelx, nely, volfrac, penal, rmin) %#ok<INUSD>
%% Define the parameters for the topology optimization problem
% If you want to change the parameters, you can do it here or pass them as arguments to the function and comment out the next line.
% There IS A way to deal with default parameters in MATLAB, but I don't want to do it(It's too long to implement and definitely not that I'm too lazy)
% So I just use the default parameters here.
nelx = 4 * 80; nely = 4 * 40; volfrac = 0.5; penal = exp(1); rmin = 3;
%% Initialize the design variable matrix x with the volume fraction
x(1:nely, 1:nelx) = volfrac;
loop = 0;
change = 1.;
%% Pre-calculate fixed and free degrees of freedom and element degrees of freedom (edof) for the finite element method
% To change boundary conditions, modify the fixeddofs and automatically set freedofs variables using setdiff.
% e.g. below is a simple fixed boundary condition where the left edge is fixed and the right edge is loaded.
fixeddofs = 1:2 * (nely + 1);
alldofs = 1:2 * (nely + 1) * (nelx + 1);
freedofs = setdiff(alldofs, fixeddofs);
edof = zeros(nely * nelx, 8);
idx = 1;

for ely = 1:nely
    for elx = 1:nelx
        n1 = (nely + 1) * (elx - 1) + ely;
        n2 = (nely + 1) * elx + ely;
        edof(idx, :) = [2 * n1 - 1, 2 * n1, 2 * n2 - 1, 2 * n2, 2 * n2 + 1, 2 * n2 + 2, 2 * n1 + 1, 2 * n1 + 2];
        idx = idx + 1;
    end
end

%% Define the material properties and stiffness matrix for the finite element analysis
% Using a simple isotropic material model with Young's modulus E and Poisson's ratio nu
% This matrix is the base stiffness matrix for a 2D linear elastic element, therefore it not be changed in each iteration.
% If any other material model is needed(Isotropic or AnIsotropic), it can be defined here:
E = 1.; nu = 0.3;
k = [1/2 - nu / 6 1/8 + nu / 8 -1/4 - nu / 12 -1/8 + 3 * nu / 8 ...
    -1/4 + nu / 12 -1/8 - nu / 8 nu / 6 1/8 - 3 * nu / 8];
% This is the stiffness matrix for a standard 4-node bilinear quadrilateral element (Q4 element).
KE = E / (1 - nu ^ 2) * ...
    [k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
%% Initialize the force vector F
ndof = 2 * (nely + 1) * (nelx + 1);
F1 = sparse(ndof, 1);
% To change the loading condition, modify the F1 vector accordingly.
% Also you can add some body forces or other loading conditions.
% e.g. to apply a unit load at the center of the right edge to simulate a simple loading condition - Cantilever Beam:
F1(2 * (nely + 1) * nelx + nely + 2, 1) = -1; % A downward vertical load at the midpoint of the right edge, which corresponds to the y-direction degree of freedom of the node at the midpoint of the right edge. If you want an upward vertical load, change it to 1.
% F1(2 * (nely + 1) * nelx + nely + 1, 1) = 1; % A horizontal load at the midpoint of the right edge, which corresponds to the x-direction degree of freedom of the node at the midpoint of the right edge. This is a rightward horizontal load, if you want a leftward horizontal load, change it to -1.

% If body forces are needed, it can be added here as well.
% If the body force varies in each loop(for example, some body forces depend on the pseudo-density x), it has to be updated accordingly. 
% Here we set F2 as a downward body force, which is positively correlated with the pseudo-density x.
body_force_scale = 5e-4; % Body force scaling factor, adjust according to the actual problem
F2 = assemble_body_force(nelx, nely, x, body_force_scale);

F = F1 + F2; % Superimpose different types of loads into the total force vector F

%% Main optimization loop
while change > 0.01 && loop < 120 % Loop until the change in design variables is small enough or a maximum number of iterations is reached
    %% Initialize the change in design variables and previous design variable matrix
    loop = loop + 1;
    x_old = x;
    %% Compute the finite element analysis
    [U] = FE(nelx, nely, x, penal, KE, freedofs, F);
    %% Calculate the objective function and its sensitivity with respect to the design variables
    [c, dc] = calc(nelx, nely, x, penal, KE, U, edof);
    %% Apply the filter to the sensitivity
    [dcn] = check(nelx, nely, rmin, x, dc);
    %% Use the Optimality Criteria method to update the design variables
	% This can be changed to other optimization algorithms such as MMA if you like, but the OC method is simple and efficient for this problem.
    [x] = OC(nelx, nely, x, volfrac, dcn);
    %% Update F2 as body force: reassemble in each loop according to the current pseudo-density x.
    F2 = assemble_body_force(nelx, nely, x, body_force_scale);
    F = F1 + F2; % Superimpose loads again
    %% Calculate the change in design variables
    % The change is the maximum absolute difference between the current and previous design variables
    change = max(abs(x(:) - x_old(:)));
    %% Show results
    show_result_per_iteration(nelx, nely, x, loop, c, change);
end
save_final_result(x); % Save the final result
end

%% ---- Function Definitions ----
%% Optimality Criteria Method for updating design variables
% I change nothing in this function. It is the same as the original top99.m.
% If you have any better version to update design variable matrix x such as mma, just replace it.
function [x_new] = OC(nelx, nely, x, volfrac, dc)
l1 = 0; l2 = 1e6; move = 0.2;

while (l2 - l1 > 1e-4)
    lmid = 0.5 * (l2 + l1);
    x_new = max(0.001, max(x - move, min(1., min(x + move, x .* sqrt(-dc ./ lmid)))));
    if sum(x_new(:)) - volfrac * nelx * nely > 0
        l1 = lmid;
    else
        l2 = lmid;
    end
end
end

%% Filter the sensitivity to ensure a smooth design update
function [dcn] = check(nelx, nely, rmin, x, dc)
%% Pre-compute the filter matrix H and its row sums Hs for efficiency and numerical stability
% The filter matrix H is used to average the sensitivities dc over a neighborhood defined by rmin.
% This is a common technique in topology optimization to avoid checkerboard patterns and ensure a smooth design.
% It uses a persistent variable to store the filter matrix H and its row sums Hs for efficiency.
% For fac is re-calculate in previous version and it never change in each iteration, this avoids recalculating the filter matrix in every iteration.
% There might be some other way to calculate the filter matrix more efficiently, but due to it calculating the filter matrix only once and reusing it, this method is sufficient enough for most cases.
persistent H Hs
if isempty(H)
    nelem = nelx * nely;
    nfilter = nelem * (2 * floor(rmin) + 1) ^ 2;
    iH = zeros(nfilter, 1, "int32");
    jH = zeros(nfilter, 1, "int32");
    sH = zeros(nfilter, 1);
    cc = 0;
    for i = 1:nelx
        for j = 1:nely
            row = (i - 1) * nely + j;
            for k = max(i - floor(rmin), 1):min(i + floor(rmin), nelx)
                for l = max(j - floor(rmin), 1):min(j + floor(rmin), nely)
                    col = (k - 1) * nely + l;
                    fac = rmin - sqrt((i - k) ^ 2 + (j - l) ^ 2);
                    if fac > 0
                        cc = cc + 1;
                        iH(cc) = row;
                        jH(cc) = col;
                        sH(cc) = fac;
                    end
                end
            end
        end
    end
    H = sparse(iH(1:cc), jH(1:cc), sH(1:cc), nelem, nelem);
    Hs = sum(H, 2);
end

%% Vectorize the design variable x and sensitivity dc for efficient computation
x_vec = x(:);
dc_vec = dc(:);
dcn_vec = (H * (x_vec .* dc_vec)) ./ (Hs .* x_vec);
dcn = reshape(dcn_vec, nely, nelx);
end

%% Finite Element Analysis function
function [U] = FE(nelx, nely, x, penal, KE, freedofs, F)
nelem = nelx * nely;
%% Pre-compute the global stiffness matrix indices iK and jK for efficiency
% Use a persistent variable to store the global stiffness matrix indices iK and jK for efficiency, avoiding recalculating them in every iteration.
persistent iK jK
if isempty(iK)
    iK = zeros(64 * nelem, 1, 'int32'); % Row Index
    jK = zeros(64 * nelem, 1, 'int32'); % Col Index
    k = 1;
    for elx = 1:nelx
        for ely = 1:nely
            n1 = (nely + 1) * (elx - 1) + ely;
            n2 = (nely + 1) * elx + ely;
            edof = [2 * n1 - 1; 2 * n1; 2 * n2 - 1; 2 * n2; 2 * n2 + 1; 2 * n2 + 2; 2 * n1 + 1; 2 * n1 + 2];
            for i = 1:8
                for j = 1:8
                    iK(k) = edof(i);
                    jK(k) = edof(j);
                    k = k + 1;
                end
            end
        end
    end
end
% This replaces assembling element stiffness matrix K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE; to calculate indices (edof, edof) for the global stiffness matrix.
% The indices iK and jK are used to construct the global stiffness matrix K in a sparse format.
% In addition, this avoids the need to loop over each element to construct the global stiffness matrix and recalculate in each iteration, which is more efficient for large problems.
% There might be some other way to calculate the indices iK and jK more efficiently, but this method is simple and efficient enough for most cases for it only require a one-time setup.
%% Vectorize the design variable x and apply the penalization
x_vec = x(:);
x_penal = x_vec .^ penal;
%% Vectorize the stiffness matrix KE and compute the global stiffness matrix entries
KE_vec = KE(:);
%% Pre-allocate the global stiffness matrix sK
% Method1: use kron
sK = kron(x_penal, KE_vec);
% Method2: use repmat and repelem
% sK_buffer = repmat(KE_vec, nelem, 1) .* repelem(x_penal, 64);
% Method3: use traditional sparse matrix construction
% sK = zeros(64 * nelem, 1);  % Values
% k = 1;
% for elx = 1:nelx
%     for ely = 1:nely
%         ke = x(ely, elx) ^ penal * KE;
%         for i = 1:8
%             for j = 1:8
%                 sK(k) = ke(i, j);
%                 k = k + 1;
%             end
%         end
%     end
% end
% These 3 methods above are equivalent, but Method1 is more efficient in terms of memory usage, Method3 is more intuitive but less efficient.
%% Initialize the global stiffness matrix K and displacement vector U
K = sparse(iK, jK, sK, 2 * (nelx + 1) * (nely + 1), 2 * (nelx + 1) * (nely + 1));
U = zeros(2 * (nely + 1) * (nelx + 1), 1);
%% Solve the system of equations K * U = F for the free degrees of freedom
U(freedofs, :) = K(freedofs, freedofs) \ F(freedofs, :);
% U(freedofs, :) = linear_system_solver(K(freedofs, freedofs), F(freedofs, :));
% U(fixeddofs, :) = 0;  % Apply boundary conditions (fixed degrees of freedom)
% It is obvious that the fixed degrees doesn't participate in solving, so why to set them to zero again?
end

%% Assemble downward body force which is positively correlated with the pseudo-density
function F2 = assemble_body_force(nelx, nely, x, body_force_scale)
% Here we distribute the body force to the four corner nodes of each element, and the magnitude of the body force is proportional to the pseudo-density x of the element (body_force_scale).
nodal_fy = zeros(nely + 1, nelx + 1);
elem_fy_each_node = -0.25 * body_force_scale * x; % The downward body force of each element is evenly distributed to the four corner nodes

nodal_fy(1:nely, 1:nelx) = nodal_fy(1:nely, 1:nelx) + elem_fy_each_node;
nodal_fy(2:nely + 1, 1:nelx) = nodal_fy(2:nely + 1, 1:nelx) + elem_fy_each_node;
nodal_fy(1:nely, 2:nelx + 1) = nodal_fy(1:nely, 2:nelx + 1) + elem_fy_each_node;
nodal_fy(2:nely + 1, 2:nelx + 1) = nodal_fy(2:nely + 1, 2:nelx + 1) + elem_fy_each_node;

ndof = 2 * (nely + 1) * (nelx + 1);
F2 = sparse(ndof, 1);
F2(2:2:end, 1) = nodal_fy(:);
end

%% Calculate objective function value c and sensitivity dc
function [c, dc] = calc(nelx, nely, x, penal, KE, U, edof)
dc = zeros(nely, nelx);
% The following line is the vectorized version of the sensitivity calculation
Ue_all = reshape(U(edof', 1), 8, []); % Extract the displacements for all elements by edof. In previous versions, it was done in a loop by Ue = U(edof(idx, :), 1);
UKU_all = sum((KE * Ue_all) .* Ue_all, 1); % Compute the element strain energy = U'KU, then sum over the element degrees of freedom
% x_vec = x(:);
% c = sum((x_vec .^ penal) .* UKU_all);
% Below it is a vectorized way to calculate the c and dc, which avoids loops by using vectorized operations
% dc_vec = -penal * (x_vec .^(penal - 1)) .* UKU_all;
% dc = reshape(dc_vec, nely, nelx);
% But the following method is efficient enough in terms of memory usage and performance and more intuitive.
idx = int32(1);
c = 0.; % Objective function value
% c calculated below which is the objective function value is not used in the optimization loop, but can be useful for debugging or monitoring progress.
for ely = 1:nely
    for elx = 1:nelx
        c = c + x(ely, elx) ^ penal * UKU_all(idx);
        dc(ely, elx) = -penal * x(ely, elx) ^ (penal - 1) * UKU_all(idx);
        idx = idx + 1;
    end
end
end

%% Show results per iteration
function show_result_per_iteration(nelx, nely, x, loop, c, change)
% Display the current iteration It.=loop, objective function value Obj.=c, volume fraction Vol. and change ch.=change.
% If best performance is desired, consider commenting out the display line and plotting (especially the plotting).
% It is useful for debugging and monitoring progress but significantly slows down the iteration speed.
disp([' It.: ' sprintf('%4i', loop) ' Obj.: ' sprintf('%10.4f', c) ...
    ' Vol.: ' sprintf('%6.3f', sum(x(:)) / (nelx * nely)) ...
    ' ch.: ' sprintf('%6.3f', change)]);
% Plot the current design variable layout using a grayscale colormap. 
% Larger x values result in lighter colors (closer to white), while smaller x values result in darker colors (closer to black).
% Displaying -x is to make the high-density areas more prominent.
colormap(gray); imagesc(-x); axis equal; axis tight; axis off;
end

%% Save the final result
function save_final_result(x)
% Here you can save the final design variable matrix x to a file, or output the result in other ways.
% For example, save as a MAT file:
save('final_design.mat', 'x');
end
