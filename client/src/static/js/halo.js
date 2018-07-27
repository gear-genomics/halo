import pako from 'pako'

$('#mainTab a').on('click', function(e) {
  e.preventDefault()
  $(this).tab('show')
})

const fileInput = document.getElementById('inputFile')
const resultLink = document.getElementById('link-results')
const resultInfo = document.getElementById('result-info')
const sampleSelect = document.getElementById('selectSample')
const resultContainer = document.getElementById('result-container')
const chartsContainer = document.getElementById('charts-container')

const submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', function() {
  resultLink.click()
  run()
})

sampleSelect.addEventListener('change', handleChangeSample)

var inputData

function run() {
  const file = fileInput.files[0]

  // TODO error handling
  if (file === undefined) return

  const fileReader = new FileReader()
  const isGzip = /\.gz$/.test(file.name)
  if (isGzip) {
    fileReader.readAsArrayBuffer(file)
  } else {
    fileReader.readAsText(file)
  }
  fileReader.onload = event => {
    let content = event.target.result
    if (isGzip) {
      content = pako.ungzip(content, { to: 'string' })
    }
    inputData = JSON.parse(content).data
    const samples = inputData.map(rec => rec.sample)
    sampleSelect.innerHTML = samples
      .map((sample, index) => `<option value="${index}">${sample}</option>`)
      .join('\n')
    showElement(resultContainer)
    vis(inputData[0])
  }
}

function handleChangeSample(event) {
  const index = Number.parseInt(event.target.value)
  vis(inputData[index])
}

function vis(data) {
  chartsContainer.innerHTML = ''
  for (const coverage of data.coverages) {
    const coverageContainer = document.createElement('div')
    const chartContainer = document.createElement('div')
    coverageContainer.appendChild(chartContainer)
    chartsContainer.appendChild(coverageContainer)
    const title = `${coverage.chromosome} – ${data.sample}`
    renderStrandedCoverageChart(chartContainer, coverage, title)
  }
}

function renderStrandedCoverageChart(container, data, title) {
  const trace1 = {
    x: data.positions,
    y: data.counts[0].values,
    name: data.counts[0].label,
    type: 'bar'
  }

  const trace2 = {
    x: data.positions,
    y: data.counts[1].values,
    yaxis: 'y2',
    name: data.counts[1].label,
    type: 'bar'
  }

  const traces = [trace1, trace2]

  const layout = {
    title: title || '',
    yaxis: {
      title: 'count',
      domain: [0.54, 1]
    },
    yaxis2: {
      domain: [0, 0.46],
      autorange: 'reversed'
    }
  }

  Plotly.newPlot(container, traces, layout)
}

function showElement(element) {
  element.classList.remove('d-none')
}

function hideElement(element) {
  element.classList.add('d-none')
}
