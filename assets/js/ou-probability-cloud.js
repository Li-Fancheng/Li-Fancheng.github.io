/* Exact OU transition density in y = x / sqrt(D/gamma), tau = gamma*t. */
(() => {
  const ns = "http://www.w3.org/2000/svg";
  document.querySelectorAll("[data-ou-cloud]").forEach((widget) => {
    if (widget.dataset.initialized) return;
    widget.dataset.initialized = "true";
    const svg = widget.querySelector("[data-ou-plot]");
    const slider = widget.querySelector("[data-ou-slider]");
    const left = 58, right = 682, top = 16, bottom = 258;
    const xmin = -4, xmax = 5, ymax = 1.45, initial = 2.5;
    const sx = (x) => left + (x - xmin) * (right - left) / (xmax - xmin);
    const sy = (y) => bottom - y * (bottom - top) / ymax;
    const add = (tag, attrs, text, parent = svg) => {
      const el = document.createElementNS(ns, tag);
      Object.entries(attrs).forEach(([key, value]) => el.setAttribute(key, value));
      if (text !== undefined) el.textContent = text;
      parent.appendChild(el);
      return el;
    };
    const densityPath = (tau, stationary = false) => {
      const mean = stationary ? 0 : initial * Math.exp(-tau);
      const variance = stationary ? 1 : -Math.expm1(-2 * tau);
      let path = "";
      for (let i = 0; i <= 720; i++) {
        const x = xmin + (xmax - xmin) * i / 720;
        const p = Math.exp(-((x - mean) ** 2) / (2 * variance)) / Math.sqrt(2 * Math.PI * variance);
        path += `${i ? "L" : "M"}${sx(x).toFixed(2)},${sy(p).toFixed(2)}`;
      }
      return path;
    };

    [0, 0.5, 1].forEach((y) => {
      add("line", {x1:left, x2:right, y1:sy(y), y2:sy(y), class:"ou-cloud-grid"});
      add("text", {x:left - 10, y:sy(y) + 5, "text-anchor":"end", class:"ou-cloud-axis-text"}, y.toFixed(1));
    });
    [-4, -2, 0, 2, 4].forEach((x) => {
      add("line", {x1:sx(x), x2:sx(x), y1:bottom, y2:bottom + 5, class:"ou-cloud-axis"});
      add("text", {x:sx(x), y:bottom + 25, "text-anchor":"middle", class:"ou-cloud-axis-text"}, String(x));
    });
    add("line", {x1:left, x2:right, y1:bottom, y2:bottom, class:"ou-cloud-axis"});
    add("line", {x1:left, x2:left, y1:top, y2:bottom, class:"ou-cloud-axis"});
    add("text", {x:370, y:304, "text-anchor":"middle", class:"ou-cloud-axis-text"}, "位置 x / ℓ");
    add("text", {x:17, y:140, transform:"rotate(-90 17 140)", "text-anchor":"middle", class:"ou-cloud-axis-text"}, "密度 ℓ p(x,t)");
    [0.05, 0.25, 0.75, 1.75, 4].forEach((tau) => {
      add("path", {d:densityPath(tau), class:"ou-cloud-reference-path"});
    });
    add("path", {d:densityPath(0, true), class:"ou-cloud-stationary-path"});
    const area = add("path", {class:"ou-cloud-area"});
    const curve = add("path", {class:"ou-cloud-current-path"});
    const meanLine = add("line", {y1:bottom, y2:top + 12, class:"ou-cloud-mean-line"});
    const delta = add("g", {});
    const dx = sx(initial);
    add("line", {x1:dx, x2:dx, y1:bottom, y2:sy(1.22), class:"ou-cloud-delta"}, undefined, delta);
    add("path", {d:`M${dx - 6},${sy(1.22) + 10}L${dx},${sy(1.22)}L${dx + 6},${sy(1.22) + 10}`, class:"ou-cloud-delta"}, undefined, delta);
    add("text", {x:dx + 12, y:sy(1.22) + 4, class:"ou-cloud-axis-text"}, "δ 初态", delta);
    const timeText = widget.querySelector("[data-ou-time]");
    const meanText = widget.querySelector("[data-ou-mean]");
    const varianceText = widget.querySelector("[data-ou-variance]");
    const explanation = widget.querySelector("[data-ou-explanation]");
    const render = () => {
      const tau = Number(slider.value);
      const mean = initial * Math.exp(-tau);
      const variance = -Math.expm1(-2 * tau);
      timeText.textContent = `γt = ${tau.toFixed(2)}`;
      meanText.textContent = `均值：${mean.toFixed(3)} ℓ`;
      varianceText.textContent = `方差：${variance.toFixed(3)} ℓ²`;
      meanLine.setAttribute("x1", sx(mean));
      meanLine.setAttribute("x2", sx(mean));
      meanLine.setAttribute("visibility", tau === 0 ? "hidden" : "visible");
      delta.setAttribute("visibility", tau === 0 ? "visible" : "hidden");
      const path = tau === 0 ? "" : densityPath(tau);
      curve.setAttribute("d", path);
      area.setAttribute("d", path ? `${path}L${right},${bottom}L${left},${bottom}Z` : "");
      explanation.textContent = tau === 0
        ? "初态严格集中在 2.5ℓ，箭头表示 δ 分布，其方差为零；灰线为后续时刻的参考形状。"
        : "蓝线表示当前分布，竖虚线标出均值。中心向原点返回，方差趋近 ℓ²；ℓ = √(D/γ)。";
      svg.setAttribute("aria-label", `OU 概率云：γt=${tau.toFixed(2)}，均值 ${mean.toFixed(3)}ℓ，方差 ${variance.toFixed(3)}ℓ²`);
    };
    slider.addEventListener("input", render);
    render();
  });
})();
