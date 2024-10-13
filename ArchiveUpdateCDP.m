% ArchiveUpdatePQ
function[Ax,Ay,A_CV] = ArchiveUpdateCDP(Ax, Ay, A_CV, Px, Py, P_CV, tol)

   for i=1:size(Py,1)
        flag = true;
        for j=size(Ay,1):-1:1
            if dominance_CV(Ay(j,:),A_CV(j),Py(i,:),P_CV(i),tol)
%             if all(Ay(j,:) <= Py(i,:))
                flag = false;
                break;
            end
            if dominance_CV(Py(i,:),P_CV(i),Ay(j,:),A_CV(j),tol)
%             if all(Py(i,:) <= Ay(j,:))
                Ax(j,:) = []; 
                Ay(j,:) = []; 
                A_CV(j) = [];
            end
        end
        if flag
            Ax = [Ax; Px(i,:)];
            Ay = [Ay; Py(i,:)];
            A_CV = [A_CV; P_CV(i)];
        end
    end
end



%---------------------Subfunctions----------------------------------------%

function dom = dominance_CV(a,a_CV,b,b_CV,eps)
%     eps = 1e-4;
    if a_CV <= eps && b_CV <= eps %BOTH ARE FEASIBLE
        dom = all(a <= b);
%         dom = prod(double(a<=b));
    elseif a_CV <= eps && b_CV > eps %a feasible, b unfeasible
        dom = true;
    elseif a_CV > eps && b_CV > eps %both unfeasible
        if a_CV<b_CV 
            dom = true;
        else
            dom = false;
        end
    elseif a_CV > eps && b_CV <= eps %a unfeasible, b feasible
        dom = false;
    end
end

function[xn,yn,yn_CV,xd,yd,yd_CV] = nondom_CV(x,y,y_CV,eps)
   n = length(y(:,1));
   I = ones(1,n);
   xd = [];
   yd = [];
   yd_CV = [];
   
   
   for i = 1:n
      for j = i+1:n
          dom = vdomw_CV(y(j,:),y_CV(j),y(i,:),y_CV(i),eps);
	  if dom==1
	     I(i) = 0;
	  elseif dom==-1
	     I(j) = 0;   
	  end
      end
   end
   
   ind  = 1;
   for i = 1:n
       if I(i) == 1
           xn(ind,:) = x(i,:);
           yn(ind,:) = y(i,:);
           yn_CV(ind) = y_CV(i);
           ind = ind + 1;
       else
           xd(end+1,:) = x(i,:);
	       yd(end+1,:) = y(i,:);
           yd_CV(end+1) = y_CV(i);
       end
   end 	 
end 

function[dom] =  vdomw_CV (v,v_CV,w,w_CV,eps)
    n = length(v);
%     eps = 1e-4;

    %CHECK IF v=w
    t = 0;  
    for i = 1:n
        if v(i) ~= w(i)  
           t = 1;
           break;
        end
    end     
    if  t==0
        dom = 0;  
        return;
    end  

    if v_CV <= eps && w_CV <= eps %BOTH ARE FEASIBLE
        t=0;
        for i = 1:n
            if v(i) < w(i) 
                t=1;
                break;
            end  
        end  
        if t==0 
            dom = -1;
            return;
        end     
        
        t=0;
        for i=1:n
            if v(i) > w(i) 
                t=1;
                break;
            end  
        end  
        if  t==0
            dom = 1;  
            return;
        end
            
        dom = 0;
    elseif v_CV <= eps && w_CV > eps %a feasible, b unfeasible
        dom = true;
    elseif v_CV > eps && w_CV > eps %both unfeasible
        if v_CV<w_CV 
            dom = true;
        else
            dom = false;
        end
    elseif v_CV > eps && w_CV <= eps %a unfeasible, b feasible
        dom = false;
    end
  
end

%-------------------------------------------------------------------------%