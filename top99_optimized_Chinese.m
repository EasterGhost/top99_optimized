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

%% 写在前面
% 对于下载此份代码的同学，非常欢迎你学习和使用它！但请务必遵守 MIT 许可证的规定，保留版权声明和许可信息，并且不要删除或修改这些信息。
% 我希望你能够规范地下载、使用和传播这份代码，尊重原作者的劳动成果，同时也欢迎你对代码进行改进和优化，并且在遵守许可证的前提下分享你的改进版本。
% 受限于此代码是对原始 top99.m 的改写，因此我并未遵循代码规范的最佳实践来编写它，例如函数划分、变量命名、注释风格等方面可能存在一些不规范的地方，但我尽量保持了代码的清晰和可读性。
% 鉴于大部分读者可能是初学者，我推荐你在以后的代码实践中逐步学习和应用更规范的代码风格和最佳实践，例如：
% · 使用更具描述性的变量名
% · 合理划分函数
% · 添加更多注释
% · 遵循一致的命名约定
% · 使用多文件结构
% · 使用版本控制系统
% · 编写单元测试
% · 规范化文档等
% 这些都是提升代码质量和可维护性的关键因素。
% 如果你有任何关于代码的疑问、建议或者改进意见，欢迎随时提交 issue，我会尽力回复和帮助你。同时也欢迎你在 GitHub 上提交 pull request 来贡献你的改进版本。
% 最后，感谢你对这份代码的关注和支持！如果你喜欢这份代码，欢迎给它点个 star！希望它能够帮助你更好地理解和应用结构拓扑优化的相关知识和技术。

%% ---- 基于最优性准则（OC）法的拓扑优化 ----
% 本代码是对原始 top99.m 的改写，是结构拓扑优化问题的 MATLAB 实现。
% 优化版采用向量化运算和持久变量以提升结构拓扑优化的效率。
% 通过索引访问和预计算矩阵（如 KE 组装全局刚度矩阵 K），提升了性能。
% 代码的核心优化思想与 top88.m 类似：利用向量化和持久变量提升多核性能并避免重复计算。持久变量不会在每次函数调用时被清除，从而节省了重复计算的时间开销。
% 但在中小规模问题（<10万单元数，更大规模未测试）上更易理解且更高效。
% 顺带一提，我也不太清楚 top88.m 为什么会稍慢一些。

%% ---- 参数说明 ----
% nelx：x 方向单元数
% nely：y 方向单元数
% volfrac：体积分数
% penal：设计变量的惩罚因子
% rmin：灵敏度滤波半径
% nelem：网格总单元数，nelx * nely
% x：设计变量矩阵，用体积分数初始化，表示各单元的伪密度
% x_old：上一轮设计变量矩阵，用于计算变化量
% x_new：使用 OC 方法更新后的设计变量矩阵
% F：力向量，这里在右边界中点施加单位竖向载荷（悬臂梁示例）
% KE：单元刚度矩阵（二维线弹性单元）
% U：位移向量
% K：全局刚度矩阵（稀疏矩阵）
% fixeddofs：固定自由度（示例为左边界固支）
% freedofs：自由自由度，为全部自由度扣除固定自由度后的集合
% alldofs：全部自由度编号集合
% edof：每个单元对应的自由度编号（单元自由度）
% c：目标函数值
% dc：目标函数对设计变量的灵敏度
% dc_vec：向量化后的灵敏度，值为-penal * (x_vec .^(penal - 1)) .* UKU_all
% loop：迭代计数器
% change：设计变量变化量，用于判断是否停止迭代
% H：滤波矩阵，用 rmin 邻域对灵敏度做加权平均
% Hs：滤波矩阵各行和，用于数值稳定性
% iK, jK：全局刚度矩阵稀疏索引，预先计算以提效
% iH, jH, sH：滤波矩阵的行列索引与权值，预先计算以提效
% idx（在 FE 中用于构造 K 的 iK/jK）：单元自由度索引，从 1 开始
% idx（在灵敏度滤波中用于构造 H 的 iH/jH）：滤波索引，从 1 开始
% Ue_all：所有单元的位移（通过 edof 提取）
% UKU_all：单元应变能，对所有单元计算 U'KU，通过 edof 提取各单元位移后批量计算
% dcn：滤波后的灵敏度
% dcn_vec：向量化的滤波灵敏度，值为(H * (x_vec .* dc_vec)) ./ (Hs .* x_vec)
% x_vec：x(:) 的向量化形式
% x_penal：惩罚后的设计变量 x .^ penal
% KE_vec：KE(:) 的向量化形式
% sK：序列化的K，用于构造稀疏全局刚度矩阵的值向量，由 KE(:) 与 x_penal 组合得到
% nfilter：滤波项总数，值为nelem * (2*floor(rmin) + 1)^2

%% ---- 运行拓扑优化问题的主函数 ----
function top99_optimized_Chinese(nelx, nely, volfrac, penal, rmin) %#ok<INUSD>
%% 定义拓扑优化问题的参数
% 你可以在这里修改参数，或以参数形式传入并注释掉下面这一行。
% MATLAB 中也可以写默认参数逻辑，但实现较繁琐（使用varargin即可，并非是我偷懒），此处直接给出默认值。
nelx = 4 * 80; nely = 4 * 40; volfrac = 0.5; penal = exp(1); rmin = 3;
%% 以体积分数初始化设计变量矩阵 x
x(1:nely, 1:nelx) = volfrac;
%% 定义一些迭代控制变量
loop = 0;
change = 1.;
%% 预计算有限元分析所需的固定/自由自由度与单元自由度（edof）
% 如需修改边界条件，请调整 fixeddofs，并通过 setdiff 自动得到 freedofs。
% 例如：左边界固定、右边界加载的简单工况（悬臂梁）。
fixeddofs = 1:2 * (nely + 1); % 左边界固定，编号为 1, 2, ..., 2*(nely+1)，每个节点的两个自由度（x 和 y 方向）均固定
alldofs = 1:2 * (nely + 1) * (nelx + 1); % 全部自由度编号集合，每个节点有两个自由度（x 和 y 方向），因此总自由度数为 2 * (nely + 1) * (nelx + 1)
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

%% 定义材料参数与单元刚度矩阵（用于有限元）
% 使用简单各向同性材料模型：杨氏模量 E、泊松比 nu。
% 下述矩阵为二维线弹性 4 节点单元的基刚度矩阵，迭代过程中不变。
% 若需其他材料模型（各向同性/各向异性），可在此处定义：
E = 1.; nu = 0.3;
k = [1/2 - nu / 6 1/8 + nu / 8 -1/4 - nu / 12 -1/8 + 3 * nu / 8 ...
	-1/4 + nu / 12 -1/8 - nu / 8 nu / 6 1/8 - 3 * nu / 8];
% 这是标准四节点双线性四边形单元（Q4）的刚度矩阵。
KE = E / (1 - nu ^ 2) * ...
	[k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
	k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
	k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
	k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
	k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
	k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
	k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
	k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
%% 初始化力向量 F
ndof = 2 * (nely + 1) * (nelx + 1);
% 如需调整载荷工况，按需修改 F。
% 也可添加体力或其他载荷。
% 例如：在右边界中点施加单位竖向载荷，以模拟悬臂梁：
F1 = sparse(ndof, 1);
F1(2 * (nely + 1) * nelx + nely + 2, 1) = -1; % 右边界中点的竖向载荷，编号为 2*(nely+1)*nelx + nely + 2，对应于右边界中点节点的 y 方向自由度，这是向下的竖向载荷，如果你想要向上的竖向载荷，请将其改为 1。
% F1(2 * (nely + 1) * nelx + nely + 1, 1) = 1; % 右边界中点的水平载荷，编号为 2*(nely+1)*nelx + nely + 1，对应于右边界中点节点的 x 方向自由度，这是向右的水平载荷，如果你想要向左的水平载荷，请将其改为 -1。

% 如果需要体力，也可在此添加。
% 若体力随迭代变化（例如依赖伪密度 x），请在循环内同步更新 F。
% 这里将 F2 设置为向下体力，且与伪密度 x 正相关：x 越大，体力越大。
body_force_scale = 5e-4; % 体力缩放因子，根据实际问题调整
F2 = assemble_body_force(nelx, nely, x, body_force_scale);

F = F1 + F2; % 将不同类型的载荷叠加到总力向量 F 中

%% 拓扑优化主循环
while change > 0.01 && loop < 120 % && loop < 50 % 当设计变量变化足够小或达到最大迭代次数时停止
	%% 初始化本轮loop值与存储上一轮设计变量
	loop = loop + 1;
	x_old = x;
	%% 有限元求解
	[U] = FE(nelx, nely, x, penal, KE, freedofs, F);
	%% 计算目标函数及其对设计变量的灵敏度
	[c, dc] = calc(nelx, nely, x, penal, KE, U, edof);
	%% 对灵敏度进行滤波
	[dcn] = check(nelx, nely, rmin, x, dc);
	%% 使用 OC 方法更新设计变量
	% 如果你喜欢，也可以改为其他优化算法，如 MMA，但 OC 方法对于这个问题来说简单且高效。
	[x] = OC(nelx, nely, x, volfrac, dcn);
	%% 更新 F2 为体力：每轮根据当前伪密度 x 重新组装。
	F2 = assemble_body_force(nelx, nely, x, body_force_scale);
	F = F1 + F2; % 重新叠加载荷
	%% 计算设计变量变化量
	% 变化量为新旧设计变量的最大绝对差
	change = max(abs(x(:) - x_old(:))); % 可以改为其他衡量方式，如 2-范数等
	%% 可视化与输出
	show_result_per_iteration(nelx, nely, x, loop, c, change);
end
save_final_result(x); % 保存最终结果
end

%% ---- 子函数定义 ----
%% 使用最优性准则（OC）更新设计变量
% 此函数与原版 top99.m 保持一致，未做改动。
% 如果你有更好的更新方法（如 MMA），可在此替换。
function [x_new] = OC(nelx, nely, x, volfrac, dc)
l1 = 0; l2 = 1e6; move = 0.2;

while (l2 - l1 > 1e-4)
	lmid = 0.5 * (l2 + l1);
	x_new = max(0.001, max(x - move, min(1., min(x + move, x .* sqrt(-dc ./ lmid))))); % OC 更新公式，经典的二分法
	if sum(x_new(:)) - volfrac * nelx * nely > 0
		l1 = lmid;
	else
		l2 = lmid;
	end
end
end

%% 灵敏度滤波，抑制棋盘格，保证设计更新的平滑
% 如果希望采用密度滤波或者其他滤波方法，可以在此处替换。但是考虑到滤波内容不同作用的结果变量不同，应当同步修改主函数中的调用与传参。
function [dcn] = check(nelx, nely, rmin, x, dc)
%% 预计算滤波矩阵 H 及其行和 Hs，以提高效率与数值稳定性
% 滤波矩阵 H 用于在 rmin 定义的邻域内对 dc 做加权平均。
% 这是拓扑优化中常用的反棋盘格与平滑手段。
% 通过持久变量缓存 H 与 Hs，避免每次迭代重建。
% 旧版本会重复计算一些不变的量，这里只构建一次并复用，已能满足多数场景的效率需求。
% 如果我们把灵敏度矩阵看作一个向量，那么新的滤波灵敏度实际上是旧灵敏度的线性组合。
% 因此，dcn_vec = A * dc_vec，其中 A = H * diag(x_vec) / (Hs .* x_vec) 是线性变换矩阵。（也许是对的？我没有仔细检查数学证明，如果我错了请纠正。但代码是正确的。）
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

%% 将 x 与 dc 向量化以便高效计算
x_vec = x(:);
dc_vec = dc(:);
dcn_vec = (H * (x_vec .* dc_vec)) ./ (Hs .* x_vec);
dcn = reshape(dcn_vec, nely, nelx);
end

%% 有限元分析（FE）
function [U] = FE(nelx, nely, x, penal, KE, freedofs, F)
nelem = nelx * nely;
%% 预先计算全局刚度矩阵的索引 iK 和 jK 以提升效率
% 显然，每个单元的自由度索引 edof 是不变的，因此全局刚度矩阵 K 的索引 iK 和 jK 也是不变的。
% 因此，通过持久变量缓存 iK 和 jK，避免每次迭代重复计算，提高效率。
persistent iK jK
if isempty(iK)
	iK = zeros(64 * nelem, 1, 'int32'); % 行索引
	jK = zeros(64 * nelem, 1, 'int32'); % 列索引
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
% 这部分替换了原 top99.m 中组装刚度矩阵的部分 K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE；计算单元刚度矩阵在全局刚度矩阵中的索引 (edof, edof)。
% 索引 iK 和 jK 变量用于构建稀疏格式的全局刚度矩阵 K。
% 此外，这避免了在每次迭代中循环遍历每个单元来构建全局刚度矩阵并重新计算，对于大规模问题更高效。
% 可能还有其他更高效计算索引 iK 和 jK 的方法，但这种方法对于大多数情况来说足够简单且高效，因为它只需要一次性设置。
% 当然你可以将其置于主函数初始化阶段，形式上不参与每次迭代，并与自由度等其他预计算量放在一起，可能能够提升一点初始化的速度。这将影响第一轮迭代的速度，但对整体影响不大。
% 这里为了代码结构清晰，仍然放在 FE 函数内。
%% 向量化设计变量 x 并计算惩罚后的设计变量 x_penal
x_vec = x(:);
x_penal = x_vec .^ penal;
%% 向量化刚度矩阵 KE
KE_vec = KE(:);
%% 计算全局刚度矩阵的内容 sK
% 方法1: 使用 Kronecker 积
sK = kron(x_penal, KE_vec);
% 方法2: 使用 repmat 和 repelem
% sK_buffer = repmat(KE_vec, nelem, 1) .* repelem(x_penal, 64); % 我不大喜欢这种写法，它并不如方法1的 Kronecker 积直观和高效。
% 方法3: 使用传统稀疏矩阵构造
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
% iK 和 jK 不变，这是因为全局刚度矩阵 K 的每个元素位置是固定的，可以直接使用这些索引访问。
% 只需更新 sK 即可反映当前设计变量 x 的影响，从而高效地构建 K。
% Kronecker积计算sK的方法通常更简洁且高效，尤其是在处理大规模问题时。相关知识可参考：
% https://en.wikipedia.org/wiki/Kronecker_product
% 方法1和2仅计算sK而不计算iK和jK，iK和jK已被预先计算并存储为持久变量（persistent）。
%% 构造全局刚度矩阵 K 与位移向量 U
K = sparse(iK, jK, sK, 2 * (nelx + 1) * (nely + 1), 2 * (nelx + 1) * (nely + 1));
% 到这里我们已经高效地构建了全局刚度矩阵 K，这完整地替换了方法3，即原99行中组装刚度矩阵的部分。
U = zeros(2 * (nely + 1) * (nelx + 1), 1); % 初始化位移向量 U
%% 在自由自由度上求解线性方程组 K * U = F
U(freedofs, :) = K(freedofs, freedofs) \ F(freedofs, :);
% U(freedofs, :) = linear_system_solver(K(freedofs, freedofs), F(freedofs, :)); % 可替换为自定义线性方程组求解器，例如GPU加速求解器
% matlab 自带的反斜杠运算符已经非常高效，它对于这样的显然正定对称实矩阵 K 采用 Cholesky 分解计算其逆矩阵，通常无需替换。
% 88行中强行规定K为对称矩阵可能并不会提高效率，因为在线弹性保守系统中，K 的计算方式是 $K^(T)=\int_V (B^(T)DB)^(T)dV = \int_V B^(T)D^(T)BdV=K$ 一定是对称的。
% 在耗散系统或非线弹性情况下（如多场耦合、非保守载荷如流体压力、具有科氏力等非保守力、非互易性变形或损伤演化等情况）K才可能不对称。
% 然而在这些情况下自然也不能强行规定K是对称的。
% 因此88行的做法并不高明且可能误导用户。
% U(fixeddofs, :) = 0;  % 施加边界条件（固定自由度）
% 固定自由度不参与求解，无需再次置零。
end

%% 组装与伪密度正相关的向下体力
function F2 = assemble_body_force(nelx, nely, x, body_force_scale)
% 这里我们将体力分布在每个单元的四个角点上，且体力大小与单元的伪密度 x 成正比(body_force_scale)。
nodal_fy = zeros(nely + 1, nelx + 1);
elem_fy_each_node = -0.25 * body_force_scale * x; % 每个单元的向下体力平均分配到四个角点

nodal_fy(1:nely, 1:nelx) = nodal_fy(1:nely, 1:nelx) + elem_fy_each_node;
nodal_fy(2:nely + 1, 1:nelx) = nodal_fy(2:nely + 1, 1:nelx) + elem_fy_each_node;
nodal_fy(1:nely, 2:nelx + 1) = nodal_fy(1:nely, 2:nelx + 1) + elem_fy_each_node;
nodal_fy(2:nely + 1, 2:nelx + 1) = nodal_fy(2:nely + 1, 2:nelx + 1) + elem_fy_each_node;

ndof = 2 * (nely + 1) * (nelx + 1);
F2 = sparse(ndof, 1);
F2(2:2:end, 1) = nodal_fy(:);
end

%% 计算目标函数值 c 和灵敏度 dc
function [c, dc] = calc(nelx, nely, x, penal, KE, U, edof)
dc = zeros(nely, nelx);
% 下行是灵敏度计算的向量化写法
Ue_all = reshape(U(edof', 1), 8, []); % 通过 edof 提取所有单元的位移；旧版本为逐单元循环 Ue = U(edof(idx, :), 1)
UKU_all = sum((KE * Ue_all) .* Ue_all, 1); % 计算单元应变能 U'KU，通过 edof 提取各单元位移后批量计算
% x_vec = x(:);
% c = sum((x_vec .^ penal) .* UKU_all);
% 下述是利用向量化避免循环来计算 c 与 dc 的做法。这段代码有点小问题，尚未修改
% dc_vec = -penal * (x_vec .^(penal - 1)) .* UKU_all;
% dc = reshape(dc_vec, nely, nelx);
% 不过当前实现已足够高效且更直观。
idx = int32(1);
c = 0.; % 目标函数值
% 以下计算 c 与 dc虽然是逐单元循环的，但是效率影响不大，因为大部分计算量已在上面向量化完成，且 matlab 可能对这部分自动地进行并行优化
% 下方计算的 c 虽在优化中未直接使用，但便于调试与监控。
for ely = 1:nely
	for elx = 1:nelx
		c = c + x(ely, elx) ^ penal * UKU_all(idx);
		dc(ely, elx) = -penal * x(ely, elx) ^ (penal - 1) * UKU_all(idx);
		idx = idx + 1;
	end
end
end

%% 展示每轮迭代结果
function show_result_per_iteration(nelx, nely, x, loop, c, change)
% 显示当前迭代 It.=loop、目标值 Obj.=c、体积分数 Vol. 与变化量 ch.=change
% 如追求速度，请注释掉结果输出与绘图（尤其是绘图）；它们有助于调试但会明显拖慢迭代（尤其是绘图）。
disp([' It.: ' sprintf('%4i', loop) ' Obj.: ' sprintf('%10.4f', c) ...
	' Vol.: ' sprintf('%6.3f', sum(x(:)) / (nelx * nely)) ...
	' ch.: ' sprintf('%6.3f', change)]);
% 绘制当前设计变量分布图，使用灰度 colormap，x 越大颜色越浅（接近白色），x 越小颜色越深（接近黑色）。
% 显示为 -x 是为了让高密度区域更明显。
colormap(gray); imagesc(-x); axis equal; axis tight; axis off;
end

%% 保存最终结果
function save_final_result(x)
% 这里可以保存最终的设计变量矩阵 x 到文件，或者以其他方式输出结果。
% 例如，保存为 MAT 文件：
save('final_design.mat', 'x');
end