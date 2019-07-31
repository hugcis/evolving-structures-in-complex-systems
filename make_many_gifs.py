"""
Script for generating many GIF images from a list of rules.
"""
import sys
import subprocess

STATES = sys.argv[1]
print(STATES)
RULES = sys.argv[2:]

for rule in RULES:
    subprocess.run(["sh", "generate_frames.sh",
                    "-n", "{}".format(STATES),
                    "{}".format(rule)])
    subprocess.run(["mv", "rule_gif/temp.gif",
                    "rule_gif/many/{}.gif".format(rule)])
