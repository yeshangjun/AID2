# AID程序
## 屏蔽的参数
1. 分析->ASCDM
2. 设置->Scale A/C Size
   ```matlab
    %Calculation Settings
    % calc = uimenu(list,'Label','计算');
    % opt(8) = uimenu(calc,'Label','Trim Mode','Checked','off',...
    %     'UserData',[0,30]);
    % opt(9) = uimenu(calc,'Label','Estimate Slipstream','Checked','off',...
    %     'UserData',[0.5,0.9]);
    % opt(10) = uimenu(calc,'Label','Multhopp''s Method','Checked','on');
   ```
3. 设置
   ```matlab
    units = uimenu(list,'Label','Units');
    % opt(14) = uimenu(units,'Label','in-oz-ft/s','Checked','off');
    % opt(17) = uimenu(units,'Label','ft-lb-kts','Checked','on');

    % opt(11) = uimenu(list,'Label','Inputs/Outputs','Checked','off');
    % opt(15) = uimenu(list,'Label','Error Check','Checked','on',...
    %     'UserData',[999,89,0]);
    % opt(16) = uimenu(list,'Label','Scroll Sensitivity','Checked','off',...
    %     'UserData',0.1);
   ```
   
