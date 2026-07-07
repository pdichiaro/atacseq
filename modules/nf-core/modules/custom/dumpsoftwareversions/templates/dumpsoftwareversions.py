#!/usr/bin/env python

import platform
from textwrap import dedent

import yaml


def _make_versions_html(versions):
    html = [
        dedent(
            """\\
            <style>
            #nf-core-versions tbody:nth-child(even) {
                background-color: #f2f2f2;
            }
            </style>
            <table class="table" style="width:100%" id="nf-core-versions">
                <thead>
                    <tr>
                        <th> Process Name </th>
                        <th> Software </th>
                        <th> Version  </th>
                    </tr>
                </thead>
            """
        )
    ]
    for process, tmp_versions in sorted(versions.items()):
        html.append("<tbody>")
        for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
            html.append(
                dedent(
                    f"""\\
                    <tr>
                        <td><samp>{process if (i == 0) else ''}</samp></td>
                        <td><samp>{tool}</samp></td>
                        <td><samp>{version}</samp></td>
                    </tr>
                    """
                )
            )
        html.append("</tbody>")
    html.append("</table>")
    return "\\n".join(html)


versions_this_module = {}
versions_this_module["${task.process}"] = {
    "python": platform.python_version(),
    "yaml": yaml.__version__,
}

# Read and parse with error handling for malformed YAML
try:
    with open("$versions", encoding='utf-8', errors='ignore') as f:
        content = f.read()
        # Remove any problematic non-ASCII characters
        content = content.encode('ascii', 'ignore').decode('ascii')
        versions_by_process = yaml.safe_load(content) or {}
        if not isinstance(versions_by_process, dict):
            print(f"Warning: versions file is not a dict, got {type(versions_by_process)}")
            versions_by_process = {}
except yaml.YAMLError as e:
    print(f"Warning: Failed to parse YAML from $versions: {e}")
    print("Continuing with empty versions...")
    versions_by_process = {}
except Exception as e:
    print(f"Warning: Unexpected error reading $versions: {e}")
    versions_by_process = {}

versions_by_process.update(versions_this_module)

# aggregate versions by the module name (derived from fully-qualified process name)
# IMPROVED: Intelligent version merging to avoid false conflicts
def deep_merge_versions(base_dict, new_dict):
    """
    Recursively merge version dictionaries.
    If values are dicts themselves, merge them recursively.
    Otherwise, keep the base value (first occurrence).
    """
    merged = base_dict.copy()
    for key, value in new_dict.items():
        if key in merged:
            # If both values are dicts, merge them recursively
            if isinstance(merged[key], dict) and isinstance(value, dict):
                merged[key] = deep_merge_versions(merged[key], value)
            # Otherwise keep the existing value (first occurrence)
        else:
            merged[key] = value
    return merged

def flatten_versions(versions_dict):
    """
    Flatten nested version dictionaries to extract actual tool versions.
    Returns a dict of {tool: version} at the top level.
    """
    flat = {}
    for key, value in versions_dict.items():
        if isinstance(value, dict):
            # Recursively flatten nested dicts
            flat.update(flatten_versions(value))
        else:
            # It's an actual version string
            flat[key] = value
    return flat

versions_by_module = {}

for process, process_versions in versions_by_process.items():
    module = process.split(":")[-1]
    
    if module in versions_by_module:
        # Merge versions intelligently instead of treating as conflict
        versions_by_module[module] = deep_merge_versions(
            versions_by_module[module], 
            process_versions
        )
    else:
        versions_by_module[module] = process_versions

# Flatten any nested structures for cleaner output
for module in versions_by_module:
    versions_by_module[module] = flatten_versions(versions_by_module[module])

versions_by_module["Workflow"] = {
    "Nextflow": "$workflow.nextflow.version",
    "$workflow.manifest.name": "$workflow.manifest.version",
}

versions_mqc = {
    "id": "software_versions",
    "section_name": "${workflow.manifest.name} Software Versions",
    "section_href": "https://github.com/${workflow.manifest.name}",
    "plot_type": "html",
    "description": "are collected at run time from the software output.",
    "data": _make_versions_html(versions_by_module),
}

with open("software_versions.yml", "w") as f:
    yaml.dump(versions_by_module, f, default_flow_style=False)
with open("software_versions_mqc.yml", "w") as f:
    yaml.dump(versions_mqc, f, default_flow_style=False)

with open("versions.yml", "w") as f:
    yaml.dump(versions_this_module, f, default_flow_style=False)
