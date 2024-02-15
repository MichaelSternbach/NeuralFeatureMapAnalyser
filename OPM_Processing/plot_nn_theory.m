function plot_nn_theory(type, distance, plot_name)
% type:      d  , d++  , d+-
% :         measured distances in units of Lambda
% plot_name: name for the exported eps file 

if nargin <2
    distance = 0:0.01:1;
end
if nargin <3
    plot_name = 'theory';
end

% common design line
switch type
    case 'd'
        n = 1.2;
        xo = 0.48;
        b = 0.047;
    case 'd++'
        n = 4.5;
        xo = 0.5;
        b = 0.05;
    case 'd+-'
        n = 1.2;
        xo = 0.46;
        b = 0.08;
    otherwise 
        disp('Unknown type')
        return
end
func = @(x) x.^n./(1+exp((x-xo)/b));
plot(distance,func(distance)/max(func(distance)),'--k','Linewidth',2,'DisplayName',plot_name)

end