/**
 * @file Favorites.cpp
 * @brief Implementation of favorite locations and their QSettings storage.
 *
 * Favorites are written as a QSettings array. Every value is stored in a readable form so
 * the settings can be inspected by hand: coordinates and the Julia constant are decimal text
 * (a double would lose the deep-zoom digits), and the set type and the Auto iteration limit
 * are words.
 */

#include "pch.h"
#include <algorithm>
#include <cmath>
#include <QSettings>
#include "Favorites.h"

using namespace Qt::StringLiterals;

namespace
{

constexpr auto favoritesKey = "favorites";           ///< Settings array holding the favorites.
constexpr auto favoritesSizeKey = "favorites/size";  ///< Length of the array; written even for an empty list.
constexpr auto descriptionKey = "description";       ///< Favorite::description.
constexpr auto centerXKey = "centerX";               ///< Favorite::centerX, as decimal text.
constexpr auto centerYKey = "centerY";               ///< Favorite::centerY, as decimal text.
constexpr auto log2ZoomKey = "log2Zoom";             ///< Favorite::log2Zoom.
constexpr auto maxIterationsKey = "maxIterations";   ///< Favorite::maxIterations, or autoName for Auto.
constexpr auto setTypeKey = "setType";               ///< Favorite::setType, as mandelbrotName or juliaName.
constexpr auto juliaRealKey = "juliaReal";           ///< Real part of Favorite::juliaConstant, as decimal text; Julia favorites only.
constexpr auto juliaImagKey = "juliaImag";           ///< Imaginary part of Favorite::juliaConstant, as decimal text; Julia favorites only.

constexpr auto mandelbrotName = "mandelbrot"_L1;  ///< Stored value of setTypeKey for QMandelbrotWidget::stMandelbrot.
constexpr auto juliaName = "julia"_L1;            ///< Stored value of setTypeKey for QMandelbrotWidget::stJulia.
constexpr auto autoName = "auto"_L1;              ///< Stored value of maxIterationsKey for QMandelbrotWidget::auto_iterations.

/**
 * @brief Read a Julia constant part from the settings.
 *
 * Favorites saved before the constant had fp128 precision hold a double, which QSettings may
 * have written in exponent form, such as "1e-05"; that form is read through a double.
 *
 * @param value Stored value.
 * @return The value, or no value if it is missing, malformed or outside [-2, 2].
 */
[[nodiscard]] std::optional<fp128_t> ReadJuliaComponent(const QVariant& value)
{
    if (const std::optional<fp128_t> parsed = ParseJuliaComponent(value.toString())) {
        return parsed;
    }

    bool isNumber = false;
    const double legacy = value.toDouble(&isNumber);
    if (!isNumber || std::abs(legacy) > 2.0) {
        return std::nullopt;
    }

    return fp128_t(legacy);
}

/**
 * @brief Build the favorites a first run starts with: ten famous Mandelbrot and Julia locations.
 *
 * Each has a fixed iteration limit: the lowest Iterations menu preset that leaves under 0.3%
 * of the view unresolved, meaning drawn black although the point escapes after more
 * iterations. Beyond that, more iterations barely change the image and only cost render time.
 * Triple Spiral Valley gets the slider's maximum of 2500 instead: points next to the parabolic
 * root of the 1/3 bulb escape so slowly that 1.6% of the view is still unresolved there.
 *
 * The Julia constants have fp128 precision. The rabbit and the Siegel disk use their exact
 * values: the period 3 center, and the parameter whose fixed point rotates by the golden mean.
 *
 * @return The default favorites in menu order.
 */
[[nodiscard]] QVector<Favorite> DefaultFavorites()
{
    constexpr auto julia = QMandelbrotWidget::stJulia;

    return {
        {.description = "Seahorse Valley", .centerX = fp128_t("-0.745428"), .centerY = fp128_t("0.113009"), .log2Zoom = 10, .maxIterations = 1024},
        {.description = "Elephant Valley", .centerX = fp128_t("0.3"), .centerY = fp128_t("0.02"), .log2Zoom = 7, .maxIterations = 2048},
        {.description = "Triple Spiral Valley", .centerX = fp128_t("-0.0888"), .centerY = fp128_t("0.6556"), .log2Zoom = 10, .maxIterations = QMandelbrotWidget::max_iterations},
        // centered on the period 3 nucleus, the largest copy of the set on the real axis
        {.description = "Mini Mandelbrot (Period 3)", .centerX = fp128_t("-1.7548776662466927"), .log2Zoom = 6, .maxIterations = 384},
        // where the period doubling bulbs along the real axis accumulate
        {.description = "Feigenbaum Point", .centerX = fp128_t("-1.4011551890920506"), .log2Zoom = 7, .maxIterations = 1536},
        // the target of the zoom sequence in Wikipedia's Mandelbrot set article
        {.description = "Wikipedia Zoom Sequence", .centerX = fp128_t("-0.743643887037158704752191506114774"), .centerY = fp128_t("0.131825904205311970493132056385139"), .log2Zoom = 16, .maxIterations = 768},
        {.description = "Douady Rabbit Julia Set", .log2Zoom = 1, .maxIterations = 128, .setType = julia, .juliaConstant = {fp128_t("-0.122561166876653619975245551820735654"), fp128_t("0.744861766619744236593170428604392367")}},
        {.description = "Basilica Julia Set", .log2Zoom = 1, .maxIterations = 128, .setType = julia, .juliaConstant = {fp128_t("-1"), fp128_t("0")}},
        {.description = "Siegel Disk Julia Set", .log2Zoom = 1, .maxIterations = 128, .setType = julia, .juliaConstant = {fp128_t("-0.390540870218400050669762600713798486"), fp128_t("-0.58678790734696875119671464305571584")}},
        {.description = "Spiral Julia Set", .log2Zoom = 1, .maxIterations = 1024, .setType = julia, .juliaConstant = {fp128_t("-0.8"), fp128_t("0.156")}},
    };
}

}  // namespace

QString Favorite::displayName() const
{
    if (description.trimmed().isEmpty()) {
        return u"(untitled)"_s;
    }

    return description;
}

QVector<Favorite> LoadFavorites()
{
    QSettings settings;
    // a first run has no array at all, while a list the user emptied is stored with size 0
    if (!settings.contains(favoritesSizeKey)) {
        return DefaultFavorites();
    }

    QVector<Favorite> favorites;

    const int count = settings.beginReadArray(favoritesKey);
    favorites.reserve(count);
    for (int i = 0; i < count; ++i) {
        settings.setArrayIndex(i);

        const std::optional<fp128_t> centerX = ParseCoordinate(settings.value(centerXKey).toString());
        const std::optional<fp128_t> centerY = ParseCoordinate(settings.value(centerYKey).toString());
        bool zoomValid = false;
        const int log2Zoom = settings.value(log2ZoomKey).toInt(&zoomValid);
        const QString setType = settings.value(setTypeKey).toString();
        // only Julia favorites store a constant
        std::optional<fp128_t> juliaReal = fp128_t {}, juliaImag = fp128_t {};
        if (setType == juliaName) {
            juliaReal = ReadJuliaComponent(settings.value(juliaRealKey));
            juliaImag = ReadJuliaComponent(settings.value(juliaImagKey));
        }
        if (!centerX || !centerY || !zoomValid || !juliaReal || !juliaImag || (setType != mandelbrotName && setType != juliaName)) {
            qWarning("Skipping favorite %d: its stored location is incomplete or malformed", i + 1);
            continue;
        }

        Favorite favorite;
        favorite.description = settings.value(descriptionKey).toString();
        favorite.centerX = *centerX;
        favorite.centerY = *centerY;
        favorite.log2Zoom = std::clamp(log2Zoom, QMandelbrotWidget::logMinZoom, QMandelbrotWidget::logMaxZoom);
        // autoName doesn't parse as a number, and neither does the missing value of a favorite
        // saved before the limit was stored; both keep the default of Auto
        bool iterationsValid = false;
        const qlonglong fixedIterations = settings.value(maxIterationsKey).toLongLong(&iterationsValid);
        if (iterationsValid && fixedIterations != QMandelbrotWidget::auto_iterations) {
            favorite.maxIterations = std::clamp<int64_t>(fixedIterations, QMandelbrotWidget::min_fixed_iterations, QMandelbrotWidget::max_iterations);
        }
        favorite.setType = (setType == juliaName) ? QMandelbrotWidget::stJulia : QMandelbrotWidget::stMandelbrot;
        if (favorite.setType == QMandelbrotWidget::stJulia) {
            favorite.juliaConstant = {*juliaReal, *juliaImag};
        }
        favorites.append(favorite);
    }
    settings.endArray();

    return favorites;
}

bool SaveFavorites(const QVector<Favorite>& favorites)
{
    QSettings settings;

    // writing a shorter array leaves the old tail entries behind, unread but still stored
    settings.remove(favoritesKey);
    settings.beginWriteArray(favoritesKey, static_cast<int>(favorites.size()));
    for (qsizetype i = 0; i < favorites.size(); ++i) {
        const Favorite& favorite = favorites[i];
        const bool isJulia = favorite.setType == QMandelbrotWidget::stJulia;

        settings.setArrayIndex(static_cast<int>(i));
        settings.setValue(descriptionKey, favorite.description);
        settings.setValue(centerXKey, ToDecimalText(favorite.centerX));
        settings.setValue(centerYKey, ToDecimalText(favorite.centerY));
        settings.setValue(log2ZoomKey, favorite.log2Zoom);
        if (favorite.maxIterations == QMandelbrotWidget::auto_iterations) {
            settings.setValue(maxIterationsKey, autoName);
        } else {
            settings.setValue(maxIterationsKey, static_cast<int>(favorite.maxIterations));
        }
        settings.setValue(setTypeKey, isJulia ? juliaName : mandelbrotName);
        if (isJulia) {
            settings.setValue(juliaRealKey, ToDecimalText(favorite.juliaConstant.real));
            settings.setValue(juliaImagKey, ToDecimalText(favorite.juliaConstant.imag));
        }
    }
    settings.endArray();
    settings.sync();

    return settings.status() == QSettings::NoError;
}
