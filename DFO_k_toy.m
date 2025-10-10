clear
load('optmization_results.mat');
labelset = labelset(1:200000,:);
data_number = length(labelset);
k_predict = []; k_max = [];
k_predict_performance = []; k_max_performance = []; k_fix_performance = [];
point_list = [[1,2,6];[1,3,6];[1,4,6];[1,5,6];[2,3,6];[2,4,6];[2,5,6];[3,4,6];[3,5,6];[4,5,6]];
% for i = 1:data_number
%     x1 = 3; x2 = 4; x3 = 6;
%     y1 = labelset(i,x1);  y2 = labelset(i,x2); y3 = labelset(i,x3);
%     a = ((x1-x2)*(y2-y3)-(y1-y2)*(x2-x3))/((x1-x2)*(x2^2-x3^2)-(x1^2-x2^2)*(x2-x3)); 
%     b = ((y1-y2)*(x2^2-x3^2)-(x1^2-x2^2)*(y2-y3))/((x1-x2)*(x2^2-x3^2)-(x1^2-x2^2)*(x2-x3)); 
%     c = y1 - a*x1^2 - b*x1;
%     if a >=0
%         if y1>y3
%             k_precit_instance = x1;
%         else
%             k_precit_instance = x3;
%         end
%     else
%         k_precit_instance = min(max(round(-b/(2*a)),1),6);
%     end
%     k_predict = [k_predict,k_precit_instance];
%     k_max = [k_max,labelset(i,end)+1];
%     k_predict_performance = [k_predict_performance,labelset(i,k_precit_instance)];
%     k_max_performance = [k_max_performance,labelset(i,labelset(i,end)+1)];
%     k_fix_performance = [k_fix_performance,labelset(i,end-1)];
% end

point_performance = [];
for j = 1:length(point_list)
    k_predict = []; k_max = [];
    k_predict_performance = []; k_max_performance = []; k_fix_performance = [];
    for i = 1:data_number
        x1 = point_list(j,1); x2 = point_list(j,2); x3 = point_list(j,3);
        y1 = labelset(i,x1);  y2 = labelset(i,x2); y3 = labelset(i,x3);
        a = ((x1-x2)*(y2-y3)-(y1-y2)*(x2-x3))/((x1-x2)*(x2^2-x3^2)-(x1^2-x2^2)*(x2-x3)); 
        b = ((y1-y2)*(x2^2-x3^2)-(x1^2-x2^2)*(y2-y3))/((x1-x2)*(x2^2-x3^2)-(x1^2-x2^2)*(x2-x3)); 
        c = y1 - a*x1^2 - b*x1;
        if a >=0
            if y1>y3
                k_precit_instance = x1;
            else
                k_precit_instance = x3;
            end
        else
            k_precit_instance = min(max(round(-b/(2*a)),1),6);
        end
        k_predict = [k_predict,k_precit_instance];
        k_max = [k_max,labelset(i,end)+1];
        k_predict_performance = [k_predict_performance,labelset(i,k_precit_instance)];
        k_max_performance = [k_max_performance,labelset(i,labelset(i,end)+1)];
        k_fix_performance = [k_fix_performance,labelset(i,end-1)];
    end
    point_performance = [point_performance,mean(k_predict_performance)];
end