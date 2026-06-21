function lp = liveplot_update_1d(lp, mode, state)
%LIVEPLOT_UPDATE_1D 刷新一维 WKB/direct 计算过程中的监测图。
%   state 应包含当前需要显示的状态量。

if isempty(lp) || ~isfield(lp, 'enabled') || ~lp.enabled
    return;
end
if nargin < 2 || isempty(mode)
    mode = lp.mode;
end
mode = lower(mode);

try
    figure(lp.fig);
    clf(lp.fig);

    theta = real(state.theta(:).');
    x = real(state.x(:));
    n = real(state.n);
    n(~isfinite(n)) = 0;
    n = max(n, 0);
    rho = real(state.rho(:));
    rho(~isfinite(rho)) = 0;
    P = real(state.P(:).');
    P(~isfinite(P)) = 0;
    P = max(P, 0);
    epsFloor = max(1e-14, max(n(:))*1e-12 + 1e-14);

    %--- 图 1：trait marginal P(theta)
    subplot(2,2,1);
    plot(theta, P, 'LineWidth', 1.5); hold on;
    if isfield(state, 'theta_m') && ~isempty(state.theta_m)
        xline(state.theta_m, '--');
    end
    if isfield(state, 'theta_peak') && ~isempty(state.theta_peak)
        xline(state.theta_peak, '-.');
    end
    grid on; xlim([theta(1), theta(end)]);
    xlabel('\theta'); ylabel('P(\theta)');
    title('Trait marginal');
    if isfield(state, 'theta_m') && isfield(state, 'theta_peak')
        legend({'P(\theta)', '\theta_m', 'selected trait'}, 'Location', 'best');
    end

    %--- 图 2：WKB 时显示 u 和 H；direct 时显示 log-scale 的 P(theta)
    subplot(2,2,2);
    if strcmp(mode, 'wkb') && isfield(state, 'u') && ~isempty(state.u)
        yyaxis left;
        uu = real(state.u(:).'); uu(~isfinite(uu)) = NaN;
        plot(theta, uu, 'LineWidth', 1.5);
        ylabel('u(\theta)');
        yyaxis right;
        if isfield(state, 'H') && ~isempty(state.H)
            HH = real(state.H(:).'); HH(~isfinite(HH)) = NaN;
            plot(theta, HH, '--', 'LineWidth', 1.2);
        end
        ylabel('H(\theta)');
        xlabel('\theta');
        grid on; xlim([theta(1), theta(end)]);
        title('Phase and effective Hamiltonian');
    else
        semilogy(theta, max(P, epsFloor), 'LineWidth', 1.5);
        grid on; xlim([theta(1), theta(end)]);
        xlabel('\theta'); ylabel('P(\theta) (log)');
        title('Trait marginal (log scale)');
    end

    %--- 图 3：总密度 rho(x)
    subplot(2,2,3);
    plot(x, rho, 'LineWidth', 1.5); hold on;
    if isfield(state, 'Kx') && ~isempty(state.Kx)
        plot(x, state.Kx(:), '--', 'LineWidth', 1.0);
        legend({'\rho(x)', 'K(x)'}, 'Location', 'best');
    end
    grid on; xlim([x(1), x(end)]);
    xlabel('x'); ylabel('\rho(x)');
    title('Spatial total density');

    %--- 图 4：n(x,theta) 热图（用 log10 提高可读性）
    subplot(2,2,4);
    imagesc(theta, x, log10(n + epsFloor));
    set(gca, 'YDir', 'normal'); colorbar;
    xlabel('\theta'); ylabel('x');
    title('log_{10}(n(x,\theta)+tiny)');

    %--- 总标题：当前步数、时间、误差和关键诊断量
    ttl = local_title_text(state, mode);
    try
        sgtitle(ttl, 'Interpreter', 'none');
    catch
        annotation(lp.fig, 'textbox', [0 0.96 1 0.04], 'String', ttl, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    end

    drawnow limitrate;

    lp.lastStep = state.step;
    lp.lastTime = state.t;
    lp.counter = lp.counter + 1;

    if isfield(state, 'saveFrame') && state.saveFrame && ~isempty(lp.frameDir)
        fname = fullfile(lp.frameDir, sprintf('%s_step_%06d.png', mode, state.step));
        try
            saveas(lp.fig, fname);
        catch
            % 某些环境下 saveas 失败时忽略，不影响主程序。
        end
    end
catch ME
    warning('liveplot:updateFailed', ...
        '在线图像刷新失败，后续将关闭 livePlot。原因: %s', ME.message);
    lp.enabled = false;
end
end

function txt = local_title_text(state, mode)
parts = {};
parts{end+1} = sprintf('%s  live monitor', upper(mode));
parts{end+1} = sprintf('step=%d', state.step);
parts{end+1} = sprintf('t=%.4g/%.4g', state.t, state.T);
if isfield(state, 'residual') && ~isempty(state.residual)
    parts{end+1} = sprintf('res=%.3e', state.residual);
end
if isfield(state, 'theta_peak') && ~isempty(state.theta_peak)
    parts{end+1} = sprintf('theta*=%.4f', state.theta_peak);
end
if isfield(state, 'theta_m') && ~isempty(state.theta_m)
    parts{end+1} = sprintf('theta_m=%.4f', state.theta_m);
end
if isfield(state, 'sigma_theta') && ~isempty(state.sigma_theta)
    parts{end+1} = sprintf('sigma=%.3e', state.sigma_theta);
end
if isfield(state, 'gammaR') && ~isempty(state.gammaR) && isfinite(state.gammaR)
    parts{end+1} = sprintf('gammaR=%.3e', state.gammaR);
end
if isfield(state, 'rho_max') && ~isempty(state.rho_max)
    parts{end+1} = sprintf('rhoMax=%.3e', state.rho_max);
end
txt = strjoin(parts, '    |    ');
end
