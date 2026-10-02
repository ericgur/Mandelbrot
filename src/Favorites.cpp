/**
 * @file Favorites.cpp
 * @brief Implementation of favorite locations and their QSettings storage.
 *
 * Favorites are written as a QSettings array. Every value is stored in a readable form so
 * the settings can be inspected by hand: coordinates are decimal text (a double would lose
 * the deep-zoom digits), and the set type and the Auto iteration limit are words.
 */

#include "pch.h"
#include <algorithm>
#include <string>
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
constexpr auto juliaRealKey = "juliaReal";           ///< Real part of Favorite::juliaConstant; Julia favorites only.
constexpr auto juliaImagKey = "juliaImag";           ///< Imaginary part of Favorite::juliaConstant; Julia favorites only.

constexpr auto mandelbrotName = "mandelbrot"_L1;  ///< Stored value of setTypeKey for QMandelbrotWidget::stMandelbrot.
constexpr auto juliaName = "julia"_L1;            ///< Stored value of setTypeKey for QMandelbrotWidget::stJulia.
constexpr auto autoName = "auto"_L1;              ///< Stored value of maxIterationsKey for QMandelbrotWidget::auto_iterations.

/**
 * @brief Round a decimal number's text to fewer fraction digits, half up.
 * @param text Decimal number with a fraction, as fp128 formats it, for example "-0.25".
 * @param end Index of the first digit to drop; it must lie past the decimal point.
 * @return The text up to @p end, rounded, for example "-0.3" for "-0.25" and an @p end of 4.
 */
[[nodiscard]] std::string RoundDecimal(const std::string& text, size_t end)
{
    std::string result = text.substr(0, end);
    if (text[end] < '5') {
        return result;
    }

    // carry through the kept digits, stepping over the decimal point
    for (size_t i = result.size(); i-- > 0;) {
        char& c = result[i];
        if (c == '-') {
            break;
        }
        if (c == '.') {
            continue;
        }
        if (c != '9') {
            ++c;
            return result;
        }
        c = '0';
    }
    // every kept digit was a 9, so the carry adds a digit: "9.96" -> "10.0"
    result.insert(result.starts_with('-') ? 1 : 0, 1, '1');

    return result;
}

/**
 * @brief Build the favorites a first run starts with: ten famous Mandelbrot and Julia locations.
 *
 * Every location was picked by rendering it with the Auto iteration limit, and they keep
 * Auto: fixed limits of 512 and up darken several of them with this palette, and Auto keeps
 * scaling when the user zooms on from there. Julia constants keep to the 10 decimals that the
 * favorites editor accepts.
 *
 * @return The default favorites in menu order.
 */
[[nodiscard]] QVector<Favorite> DefaultFavorites()
{
    constexpr auto julia = QMandelbrotWidget::stJulia;
    constexpr auto autoLimit = QMandelbrotWidget::auto_iterations;

    return {
        {.description = "Seahorse Valley", .centerX = fp128_t("-0.745428"), .centerY = fp128_t("0.113009"), .log2Zoom = 10, .maxIterations = autoLimit},
        {.description = "Elephant Valley", .centerX = fp128_t("0.3"), .centerY = fp128_t("0.02"), .log2Zoom = 7, .maxIterations = autoLimit},
        {.description = "Triple Spiral Valley", .centerX = fp128_t("-0.0888"), .centerY = fp128_t("0.6556"), .log2Zoom = 10, .maxIterations = autoLimit},
        // centered on the period 3 nucleus, the largest copy of the set on the real axis
        {.description = "Mini Mandelbrot (Period 3)", .centerX = fp128_t("-1.7548776662466927"), .log2Zoom = 6, .maxIterations = autoLimit},
        // where the period doubling bulbs along the real axis accumulate
        {.description = "Feigenbaum Point", .centerX = fp128_t("-1.4011551890920506"), .log2Zoom = 7, .maxIterations = autoLimit},
        // the target of the zoom sequence in Wikipedia's Mandelbrot set article
        {.description = "Wikipedia Zoom Sequence",
         .centerX = fp128_t("-0.743643887037158704752191506114774"),
         .centerY = fp128_t("0.131825904205311970493132056385139"),
         .log2Zoom = 16,
         .maxIterations = autoLimit},
        {.description = "Douady Rabbit Julia Set", .log2Zoom = 1, .maxIterations = autoLimit, .setType = julia, .juliaConstant = {-0.123, 0.745}},
        {.description = "Basilica Julia Set", .log2Zoom = 1, .maxIterations = autoLimit, .setType = julia, .juliaConstant = {-1.0, 0.0}},
        {.description = "Siegel Disk Julia Set", .log2Zoom = 1, .maxIterations = autoLimit, .setType = julia, .juliaConstant = {-0.3905408702, -0.5867879073}},
        {.description = "Spiral Julia Set", .log2Zoom = 1, .maxIterations = autoLimit, .setType = julia, .juliaConstant = {-0.8, 0.156}},
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

const QRegularExpression& CoordinateRegularExpression()
{
    static const QRegularExpression pattern(uR"([+-]?(\d{1,2}(\.\d*)?|\.\d+))"_s);
    return pattern;
}

std::optional<fp128_t> ParseCoordinate(const QString& text)
{
    static const QRegularExpression anchored(QRegularExpression::anchoredPattern(CoordinateRegularExpression().pattern()));
    if (!anchored.match(text).hasMatch()) {
        return std::nullopt;
    }

    return fp128_t(text.toStdString());
}

QString CoordinateToString(const fp128_t& value)
{
    // fp128 prints every meaningful digit, so a value typed as "-0.2" would come back as
    // "-0.200000000000000000000000000000000001", 0.2 having no exact binary form. The fewest
    // digits that parse back to the same value keep such numbers as they were typed, and
    // checking each candidate with the parser also makes the round trip exact where possible.
    const std::string full = static_cast<std::string>(value);
    const size_t point = full.find('.');
    if (point == std::string::npos) {
        return QString::fromStdString(full);
    }

    for (size_t end = point + 2; end < full.size(); ++end) {
        const std::string candidate = RoundDecimal(full, end);
        if (fp128_t(candidate) == value) {
            return QString::fromStdString(candidate);
        }
    }

    return QString::fromStdString(full);
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
        if (!centerX || !centerY || !zoomValid || (setType != mandelbrotName && setType != juliaName)) {
            qWarning("Skipping favorite %d: its stored location is incomplete or malformed", i + 1);
            continue;
        }

        Favorite favorite;
        favorite.description = settings.value(descriptionKey).toString();
        favorite.centerX = *centerX;
        favorite.centerY = *centerY;
        favorite.log2Zoom = std::clamp(log2Zoom, static_cast<int>(QMandelbrotWidget::logMinZoom), static_cast<int>(QMandelbrotWidget::logMaxZoom));
        // autoName doesn't parse as a number, and neither does the missing value of a favorite
        // saved before the limit was stored; both keep the default of Auto
        bool iterationsValid = false;
        const qlonglong fixedIterations = settings.value(maxIterationsKey).toLongLong(&iterationsValid);
        if (iterationsValid && fixedIterations != QMandelbrotWidget::auto_iterations) {
            favorite.maxIterations = std::clamp<int64_t>(fixedIterations, QMandelbrotWidget::min_fixed_iterations, QMandelbrotWidget::max_iterations);
        }
        favorite.setType = (setType == juliaName) ? QMandelbrotWidget::stJulia : QMandelbrotWidget::stMandelbrot;
        if (favorite.setType == QMandelbrotWidget::stJulia) {
            favorite.juliaConstant = {settings.value(juliaRealKey, favorite.juliaConstant.real()).toDouble(),
                                      settings.value(juliaImagKey, favorite.juliaConstant.imag()).toDouble()};
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
        settings.setValue(centerXKey, CoordinateToString(favorite.centerX));
        settings.setValue(centerYKey, CoordinateToString(favorite.centerY));
        settings.setValue(log2ZoomKey, favorite.log2Zoom);
        if (favorite.maxIterations == QMandelbrotWidget::auto_iterations) {
            settings.setValue(maxIterationsKey, autoName);
        } else {
            settings.setValue(maxIterationsKey, static_cast<int>(favorite.maxIterations));
        }
        settings.setValue(setTypeKey, isJulia ? juliaName : mandelbrotName);
        if (isJulia) {
            settings.setValue(juliaRealKey, favorite.juliaConstant.real());
            settings.setValue(juliaImagKey, favorite.juliaConstant.imag());
        }
    }
    settings.endArray();
    settings.sync();

    return settings.status() == QSettings::NoError;
}
