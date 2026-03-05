// pkgdown/extra.js
document.addEventListener("DOMContentLoaded", function() {

  // Function to swap the logo filename while preserving relative paths
  function updateLogos(theme) {
    const logos = document.querySelectorAll("img[src$='logo.png'], img[src$='logo_dark.png']");
    logos.forEach(img => {
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
});
