/* 
Description: interactive D3.js visualization for a lattice protein folding model.

How to run locally:
    1. cd into the repository
    2. run: python3 -m http.server 8000
    3. open: http://localhost:8000
*/

const data = {
    // len=8 chain
    nodes: [
        { id: "0", type: "H" },
        { id: "1", type: "P" },
        { id: "2", type: "H" },
        { id: "3", type: "H" },
        { id: "4", type: "H" },
        { id: "5", type: "H" },
        { id: "6", type: "H" },
        { id: "7", type: "H" },
        { id: "8", type: "P" },
        { id: "9", type: "H" },
        { id: "10", type: "H" },
        { id: "11", type: "P" },
        { id: "12", type: "H" },
        { id: "13", type: "H" },
        { id: "14", type: "P" },
        { id: "15", type: "H" }
    ],
    links: [
        { source: "0", target: "1" },
        { source: "1", target: "2" },
        { source: "2", target: "3" },
        { source: "3", target: "4" },
        { source: "4", target: "5" },
        { source: "5", target: "6" },
        { source: "6", target: "7" },
        { source: "7", target: "8" },
        { source: "8", target: "9" },
        { source: "9", target: "10" },
        { source: "10", target: "11" },
        { source: "11", target: "12" },
        { source: "12", target: "13" },
        { source: "13", target: "14" },
        { source: "14", target: "15" }
  ]
};


function renderChart(){

    // Strengths of each variable
    const LINK_STRENGTH = 2.0;
    const MOUSE_STRENGTH = 2.0;
    const CELL_STRENGTH = 0.2;
    const ALPHA = 0.5;
    
    // Specify the dimensions of the chart.
    const width = 600;
    const height = 420;

    // Specify the color scale.
    const color = d3.scaleOrdinal(d3.schemeCategory10);

    // The force simulation mutates links and nodes, so create a copy
    // so that re-evaluating this cell produces the same result.
    const links = data.links.map(d => ({...d}));
    const nodes = data.nodes.map(d => ({...d}));

    // Create a simulation with several forces.
    const simulation = d3.forceSimulation(nodes)
        .force("collide", d3.forceCollide().radius(9)) // HAS STRENGTH PARAMETER
        .force("link", d3.forceLink(links).id(d => d.id).strength(LINK_STRENGTH));

    // Create the SVG container.
    const svg = d3.select("#vis")
        .append("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("viewBox", [-width / 2, -height / 2, width, height])
        .attr("style", "max-width: 100%; height: auto; user-select: none;");


    // Make grid 
    const cellsize = 30; // distance between lines
    const gridcolor = "#555"; // dark gray gridlines
    const grid = svg.append("g")
        .attr("stroke", gridcolor)
        .attr("stroke-width", 0.5);

    // Vertical lines
    for (let x = -width / 2; x <= width / 2; x += cellsize) {
        grid.append("line")
        .attr("x1", x)
        .attr("y1", -height / 2)
        .attr("x2", x)
        .attr("y2", height / 2);
    }
    
    // Horizontal lines
    for (let y = -height / 2; y <= height / 2; y += cellsize) {
        grid.append("line")
        .attr("x1", -width / 2)
        .attr("y1", y)
        .attr("x2", width / 2)
        .attr("y2", y);
    }
    
    const maxScore = 9;   // highest score
    const boxSize = 25;   // box size (smaller than before)
    const thresholds = [4, 8, 9]; // star cutoffs
    const padding = 10;
    
    // Group for boxes + stars
    const scoreGroup = svg.append("g")
        .attr("transform", `translate(${-width/2 + 22}, ${-height/2 + 45})`);

    scoreGroup.append("rect")
        .attr("x", -padding)
        .attr("y", -padding - 25) // extra space for stars
        .attr("width", maxScore * boxSize + padding * 2)
        .attr("height", boxSize + padding * 2 + 20) // extra space for stars
        .attr("fill", "white")
        .attr("stroke", "#ccc")
        .attr("rx", 8)  // rounded corners
        .attr("ry", 8);

    
    // Draw score boxes
    const boxes = [];
    for (let i = 1; i <= maxScore; i++) {
        const box = scoreGroup.append("rect")
        .attr("x", (i - 1) * boxSize)
        .attr("y", 0)
        .attr("width", boxSize)
        .attr("height", boxSize)
        .attr("fill", "#eee")
        .attr("stroke", "#999")
        .style("pointer-events", "none"); // prevents highlight/select
    
        boxes.push(box);
    }
    
    const stars = []; // store star elements

    thresholds.forEach(t => {
        const xPos = (t - 1) * boxSize + boxSize / 2;
    
        const star = scoreGroup.append("text")
        .attr("x", (t - 1) * boxSize + boxSize) // center over that box
        .attr("y", -15)
        .attr("text-anchor", "middle")
        .attr("font-size", 20)
        .attr("fill", "grey")   // start as grey
        .attr("class", "score-star") // <- give class
        .style("pointer-events", "none")
        .text("★");
    
        stars.push(star);
    
        // Vertical line at threshold
        scoreGroup.append("line")
        .attr("x1", xPos+12.5)
        .attr("x2", xPos+12.5)
        .attr("y1", -10)
        .attr("y2", boxSize)
        .attr("stroke", "#333")
        .attr("stroke-width", 2)
        .style("pointer-events", "none");
    });
    
    // function to update boxes and stars based on HH contacts
    function updateScore() {
    const hhContacts = HHContacts();
    
    // Fill boxes green for each HH contact
    for (let i = 0; i < maxScore; i++) {
        boxes[i].attr("fill", i < hhContacts ? "forestgreen" : "#eee");
    }
    
    // Update stars based on thresholds
    thresholds.forEach((t, i) => {
        stars[i].attr("fill", hhContacts >= t ? "gold" : "grey");
    });
    }
    
    
    // Function to mark box green when contact is made
    function markContact(score) {
        if (score >= 1 && score <= maxScore) {
        boxes[score - 1].attr("fill", "forestgreen");
        }
    }

    // Add a line for each link, and a circle for each node.
    const link = svg.append("g")
        .attr("stroke", "#999")
        .attr("stroke-opacity", 0.5)
        .selectAll("line")
        .data(links)
        .join("line")
        .attr("stroke-width", 2); //d => Math.sqrt(d.value));

    const node = svg.append("g")
        .attr("stroke", "#fff")
        .attr("stroke-width", 1.5)
        .selectAll("circle")
        .data(nodes)
        .join("circle")
        .attr("r", (d, i) => i === 0 ? 8 : 5)  // Make first node radius 8, others 5
        .attr("fill", function(d) {
            if (d.type == 'H') {
                return "darkorange";
            } else {
                return "steelblue";
            }
            });
    
    node.append("title")
        .text(d => d.id);

    // Use forces to drag the entire chain by each node
    let dragnode = null;
    let mouse = null;

    nodes.forEach(function(node) { node.x = node.id * cellsize + cellsize / 2; node.y = cellsize / 2; });

    node.on("mousedown", function() { dragnode = this; d3.select(dragnode).style("stroke", "#000");});
    d3.select(window).on("mouseup", function() { d3.select(dragnode).style("stroke", "#fff"); dragnode = null; });
    svg.on("mousemove", function(event) { mouse = d3.pointer(event); });

    function forceTowardMouse(alpha) {
        return function(node) {
        if (mouse !== null && dragnode !== null && node.id === dragnode.__data__.id) {
            const targetX = mouse[0];
            const targetY = mouse[1];

            const dx = targetX - node.x;
            const dy = targetY - node.y;
            node.vx += dx * alpha * MOUSE_STRENGTH;
            node.vy += dy * alpha * MOUSE_STRENGTH;
        }
        };
    }

    function forceTowardCell(alpha) {
        return function(node) {
        const targetX = Math.floor(node.x / cellsize) * cellsize + cellsize / 2;
        const targetY = Math.floor(node.y / cellsize) * cellsize + cellsize / 2;
        
        const dx = targetX - node.x;
        const dy = targetY - node.y;
        node.vx += dx * alpha * CELL_STRENGTH;
        node.vy += dy * alpha * CELL_STRENGTH;
        };
    }
    
    // Find adjacent hydrophobic nodes
    function adjacent(a, b, tolerance = 2) {
        const dx = Math.abs(a.x - b.x);
        const dy = Math.abs(a.y - b.y);
        // Instead of exact match, allow some wiggle room (±2 px)
        return (
        (Math.abs(dx - cellsize) < tolerance && dy < tolerance) ||
        (Math.abs(dy - cellsize) < tolerance && dx < tolerance)
        );
    }

    // Calculates the number of adjacent contact of hydrophobic nodes
    function HHContacts() {
        let count = 0;
        for (let i = 0; i < nodes.length; i++) {
        const a = nodes[i];
        if (a.type !== 'H') continue;
    
        for (let j = 0; j < nodes.length; j++) {
            const b = nodes[j];
            if (j === i || b.type !== 'H') continue;
    
            // neighbors on chain
            if (Math.abs(i - j) === 1) continue;
    
            if (adjacent(a, b)) {
            count++;
            }
        }
        }
    
        return count /2; // each contact is counted twice (i-j and j-i)
    }

    // Make links red if invalid 
    function isInvalidLink(link) {
        // Check if the link is not 90° or 180°
        const dx = Math.abs(link.source.x - link.target.x);
        const dy = Math.abs(link.source.y - link.target.y);
        const isStraight = (dx < 1 && Math.abs(dy - cellsize) < 1) || (dy < 1 && Math.abs(dx - cellsize) < 1);
        if (!isStraight) return true;
    
        // Count how many nodes occupy each cell
        const cellMap = {};
        nodes.forEach(n => {
        const key = `${Math.floor(n.x / cellsize)},${Math.floor(n.y / cellsize)}`;
        cellMap[key] = (cellMap[key] || 0) + 1;
        });
    
        // Check if multiple nodes in one cell 
        const sourceKey = `${Math.floor(link.source.x / cellsize)},${Math.floor(link.source.y / cellsize)}`;
        const targetKey = `${Math.floor(link.target.x / cellsize)},${Math.floor(link.target.y / cellsize)}`;
        if (cellMap[sourceKey] > 1 || cellMap[targetKey] > 1) return true;
    
        return false;
    }
    
    // Set the position attributes of links and nodes each time the simulation ticks.
    simulation.on("tick", () => {
        simulation.alpha(ALPHA);

        simulation.nodes().forEach(forceTowardMouse(simulation.alpha()));
        simulation.nodes().forEach(forceTowardCell(simulation.alpha()));
        
        link
            .attr("x1", d => d.source.x)
            .attr("y1", d => d.source.y)
            .attr("x2", d => d.target.x)
            .attr("y2", d => d.target.y)
            .attr("stroke", d => isInvalidLink(d) ? "red" : "#999") // invalid link
            .attr("stroke-opacity", d => isInvalidLink(d) ? 0.8 : 0.3);

        node
            .attr("cx", d => d.x)
            .attr("cy", d => d.y);

        // score based on HH contacts
        const score = HHContacts();
        const currentPoints = HHContacts();

        const hhLinks = svg.selectAll(".hh-contact")
        .data(nodes.flatMap((a, i) =>
            nodes.map((b, j) => {
            if (
                i < j &&
                a.type === 'H' &&
                b.type === 'H' &&
                Math.abs(i - j) > 1 &&
                adjacent(a, b)
            ) {
                return { source: a, target: b };
            } else {
                return null;
            }
            }).filter(d => d)
        ));
        
        hhLinks.enter()
        .append("line")
        .attr("class", "hh-contact")
        .attr("stroke", "orange")
        .attr("stroke-width", 2)
        .attr("stroke-opacity", 1.5)
        .attr("stroke-dasharray", "3.5, 3.5")
        .each(function() { d3.select(this).lower(); })
        .merge(hhLinks)
        .attr("x1", d => d.source.x)
        .attr("y1", d => d.source.y)
        .attr("x2", d => d.target.x)
        .attr("y2", d => d.target.y);
        
        hhLinks.exit().remove();

        updateScore();
        
    });

    return svg.node();
}


renderChart();