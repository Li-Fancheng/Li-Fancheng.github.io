// Follow relocated sections and reveal the calculation containing an anchor.
(() => {
  const element = document.getElementById("qbio-legacy-anchors");
  if (!element) return;
  const anchors = JSON.parse(element.textContent);
  const followFragment = () => {
    let fragment;
    try {
      fragment = decodeURIComponent(window.location.hash.slice(1));
    } catch {
      return;
    }
    if (!fragment) return;
    const localTarget = document.getElementById(fragment);
    if (localTarget) {
      let revealed = false;
      for (let ancestor = localTarget.parentElement; ancestor; ancestor = ancestor.parentElement) {
        if (ancestor.tagName === "DETAILS" && !ancestor.open) {
          ancestor.open = true;
          revealed = true;
        }
      }
      if (revealed) requestAnimationFrame(() => localTarget.scrollIntoView());
      return;
    }
    if (Object.prototype.hasOwnProperty.call(anchors, fragment)) {
      const target = new URL(anchors[fragment], window.location.href);
      target.search = window.location.search;
      window.location.replace(target.href);
    }
  };
  followFragment();
  window.addEventListener("hashchange", followFragment);
})();
