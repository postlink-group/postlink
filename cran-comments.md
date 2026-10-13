## Submission Notes

This is a update of the archived postlink 0.1.1 addressing the CRAN check issues and also resolving the NOTE regarding the '-Wa,-mbig-obj' flag.

## Changes in 0.1.2

-   Updated all array declarations across `.stan` files in `inst/stan/` to use the `array[N] type var;` syntax required by `stanc3`.
-   Bumped minimum required versions of `rstan` and `StanHeaders` to `>= 2.26.0` in `DESCRIPTION`.
-   Regenerated C++ model headers via `rstantools`.
-   Removed the non-portable flag(s): '-Wa,-mbig-obj'.

## Test environments

-   GitHub Actions: macOS-latest (release)
-   GitHub Actions: windows-latest (release)
-   GitHub Actions: ubuntu-latest (devel, release, oldrel-1)
-   win-builder (R-devel) (https://win-builder.r-project.org/)
-   Local Windows 11 (64-bit), R 4.5.0: `R CMD check --as-cran`

## R CMD check results

0 errors \| 0 warnings \| 1 notes

checking CRAN incoming feasibility ... [10s] **NOTE**

Maintainer: 'Priyanjali Bukke \<postlink.group\@gmail.com\>' New submission Package was archived on CRAN Possibly misspelled words in DESCRIPTION: Gutman (22:753) PLDA (22:84) Slawski (22:669) Vo (22:843) al (22:614, 22:680, 22:763, 22:849) et (22:611, 22:677, 22:760, 22:846) unlinked (22:316) CRAN repository db overrides: X-CRAN-Comment: Archived on 2026-09-27 as issues were not corrected in time. Found the following (possibly) invalid URLs: URL: <https://hdl.handle.net/10779/uow.27788247> From: DESCRIPTION Status: 202 Message: Accepted

**Explanation:**

-   The maintainer information is accurate. I would like to confirm that the DESCRIPTION file lists my name and email address as the single designated maintainer, and that the email address is not a mailing list. I do not distribute CRAN messages, and the address is used solely for maintaining our software packages (as my university email account will expire soon). I have also updated the Gmail account profile name accordingly ("postlink maintainer"), in case that is helpful.
-   All of the CRAN check issues have been addressed (please see above for the changes and these changes are also present in the resubmissions on 09-04-2026 and 09-23-2026 before the 2026-09-27 deadline in case for reference). The package should now align with the updates to rstantools and StanHeaders.
-   Thank you for these checks. We have reviewed the words and URL noted, and they appear as intended.

## Reverse dependencies

There are currently no reverse dependencies for this package.
