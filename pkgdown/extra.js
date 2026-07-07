// pkgdown/extra.js
document.addEventListener("DOMContentLoaded", function() {

  // Identify all logo elements on the page
  const logos = document.querySelectorAll("img[src$='logo.png'], img[src$='logo_dark.png'], .navbar-brand img");

  // Check if we are on an article/vignette page
  const isArticlePage = window.location.pathname.includes("/articles/") ||
                        document.querySelector(".template-article, .article-header, body.pkgdown-article");

  if (isArticlePage) {
    // Completely remove the logo for articles

    // Delete the element from the DOM.
    logos.forEach(img => img.remove());

    // Fallback: CSS just in case the navbar loads asynchronously later
    const style = document.createElement("style");
    style.type = "text/css";
    style.innerHTML = `
      .navbar-brand img,
      .navbar-brand .logo,
      a.navbar-brand > img {
        display: none !important;
        visibility: hidden !important;
        width: 0px !important;
        height: 0px !important;
      }
    `;
    document.head.appendChild(style);

  } else {
    // Run the Dark/Light theme swapper

    function updateLogos(theme) {
      logos.forEach(img => {
        if (!img.src) return; // Prevent errors if src is missing
        if (theme === "dark") {
          img.src = img.src.replace("logo.png", "logo_dark.png");
        } else {
          img.src = img.src.replace("logo_dark.png", "logo.png");
        }
      });
    }

    // Set the correct logo immediately on page load
    const currentTheme = document.documentElement.getAttribute("data-bs-theme") || "light";
    updateLogos(currentTheme);

    // Watch the HTML tag for changes when the user clicks the Light/Dark toggle
    const observer = new MutationObserver(function(mutations) {
      mutations.forEach(function(mutation) {
        if (mutation.attributeName === "data-bs-theme") {
          const newTheme = document.documentElement.getAttribute("data-bs-theme");
          updateLogos(newTheme);
        }
      });
    });

    observer.observe(document.documentElement, { attributes: true });
  }
});
