# Phasor dynamics

This directory contains an implementation of a system model using phasor
dynamics. Aside from all of the modeling aspects, one additional part of
note is the experimental JSON parsing functionality which can take in a JSON
object (described in `INPUT_FORMAT.md`) and yields a system model. This
is implemented with the [nlohmann/json](https://github.com/nlohmann/json)
parser by adding a `from_json` function for each data structure of note
with the signature `void from_json(const nlohmann::json&, SomeType&)`.

The API reference for the `nlohmann::json` type passed in can be found
[here](https://json.nlohmann.me/api/basic_json) (as well as for the rest
of the library using the navigation on that page). Documentation for some
particular methods of note on this type are here:

- `.at`: https://json.nlohmann.me/api/basic_json/at/
- `.items`: https://json.nlohmann.me/api/basic_json/items/
- `.get_to`: https://json.nlohmann.me/api/basic_json/get_to/

One final thing to note is that the existing parsing code avoids the use
of the `operator[]` implementation on `nlohmann::json` due its lack of
bounds checking. Prefer the use of `.at`.
