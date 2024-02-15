function plot_count_sd_theory(A, plot_name)
% areas:     area in units of Lambda where the pinwheels are counted
% counts:    number of pinwheels in that area
% plot_name: name for the exported eps file 

if nargin < 1
    A = [0.3 20];
end
if nargin<2
    plot_name = 'theory';
end

% common design line
gamma = 0.4;
c = 1.05;
rho = pi;
func = @(x)c*(rho./x).^gamma;


loglog(A,func(A),'--k','Linewidth',2,'DisplayName',plot_name)

end